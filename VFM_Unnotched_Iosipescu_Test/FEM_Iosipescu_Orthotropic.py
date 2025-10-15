"""
Virtual Fields Method (VFM)
---------------------------------------------------------------
Case Study: Unnotched Iosipescu Test for Orthotropic Material

Based on:
"The Virtual Fields Method: Extracting Constitutive Mechanical
Parameters from Full-Field Deformation Measurements"
by F. Pierron and M. Grédiac

Developed by:
José Xavier
Universidade NOVA de Lisboa, NOVA FCT, UNIDEMI
https://userweb.fct.unl.pt/~jmc.xavier/

Course Context:
CISM Advanced School
"Image-Based Mechanics: An Overview of Experimental and
Numerical Approaches" — Udine, Italy, 6–10 October 2025

Coordinators: Julien Réthoré and José Xavier
Lecture: "The Virtual Fields Method — Hands-On Session"

More information:
https://cism.it/en/activities/courses/C2516/

VARIABLE NAMING CONVENTION:
---------------------------------------------------------------
Geometry:
  L      = Length of the region of interest (ROI) at specimen center
  W      = Specimen width
  Th     = Specimen thickness
  tL     = Specimen total length
  IBC    = Inner boundary condition distance
  OBC    = Outer boundary condition distance
  ang    = Off-axis angle (degrees)

Mesh discretization:
  LdH1   = Number of horizontal divisions (outer wedge zones)
  LdH2   = Number of horizontal divisions (transition zones)
  LdH3   = Number of horizontal divisions (central ROI - S2 region)
  LdV    = Number of vertical divisions

Material properties (Orthotropic):
  EL, ER, ET       = Young's moduli (L=Longitudinal, R=Radial, T=Tangential)
  GLR, GLT, GRT    = Shear moduli
  CPRT, CPTL, CPLR = Poisson's ratios (e.g., CPRT = nu_RT)

Boundary conditions:
  uy     = Prescribed vertical displacement (m)
  condF  = Boundary condition type: 1=base, 2=iterative

Solution:
  SOL_el_nod = Solution output type: 1=element centroid, 2=nodal

Regions:
  S1     = Left wedge region (element component)
  S2     = Central region of interest (element component)
  S3     = Right wedge region (element component)

Input/Output:
  Input:  inputopt.dat (23 design variables)
  Output: scalarsFE.csv (global results: P, geometry)
          FEM2VFM.csv (element solution fields)
          nodal_fieldsFE.csv (nodal solution fields)
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Dict, List, Tuple, Literal

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.patches import Polygon

warnings.filterwarnings("ignore")

# =============================================================================
# CONSTANTS
# =============================================================================
TOLERANCE = 1e-9  # Geometric tolerance for BC selection (mm)
JACOBIAN_MIN = 1e-12  # Minimum acceptable Jacobian determinant
GAUSS_COORD = 1.0 / np.sqrt(3.0)  # 2×2 Gauss integration coordinate

# =============================================================================
# DATA STRUCTURES
# =============================================================================
@dataclass
class IosipescuParams:
    """Iosipescu test specimen parameters."""

    # Geometry (mm)
    ang: float  # Material off-axis angle (degrees)
    L: float  # Length of S2 region (gauge section)
    W: float  # Specimen width
    Th: float  # Specimen thickness
    tL: float  # Total specimen length
    OBC: float  # Outer grip length (SHORT side)
    IBC: float  # Inner grip length (LONG side)

    # Material properties (MPa for elastic moduli and shear moduli)
    EL: float  # Longitudinal modulus
    ER: float  # Radial modulus
    ET: float  # Tangential modulus
    CPRT: float  # Poisson's ratio νRT
    CPTL: float  # Poisson's ratio νTL
    CPLR: float  # Poisson's ratio νLR
    GRT: float  # Shear modulus GRT
    GLT: float  # Shear modulus GLT
    GLR: float  # Shear modulus GLR

    # Mesh divisions
    LdH3: int  # Divisions in S2 region (x-direction)
    LdV: int  # Divisions in y-direction
    LdH1: int  # Divisions in outer grip regions
    LdH2: int  # Divisions in transition regions

    # Loading (mm)
    uy: float  # Prescribed vertical displacement

    def __post_init__(self):
        """Validate inputs after initialization."""
        self._validate()

    def _validate(self):
        """Validate units and physical constraints."""
        # Check units (MPa scale)
        if not (100 < self.EL < 1e6):
            warnings.warn(f"EL={self.EL} seems unusual. Expected MPa (1000-50000)?")
        if not (0 < self.uy < 100):
            warnings.warn(f"uy={self.uy} seems unusual. Expected mm (0.1-10)?")

        # Check Poisson's ratios
        if not (0 < self.CPLR < 0.5):
            raise ValueError(f"Invalid CPLR={self.CPLR}. Must be in (0, 0.5)")

        # Check geometry
        if self.OBC >= self.IBC:
            raise ValueError(f"OBC={self.OBC} must be < IBC={self.IBC}")
        if self.L <= 0 or self.W <= 0 or self.Th <= 0:
            raise ValueError("Geometry dimensions must be positive")

        print(f"✓ Parameters validated (EL={self.EL:.0f} MPa, uy={self.uy:.3f} mm)")


# =============================================================================
# MAIN FEM SOLVER CLASS
# =============================================================================
class IosipescuFEM:
    """
    Finite element solver for Iosipescu shear test on orthotropic materials.

    Features:
    - 4-node quad elements (PLANE182-like)
    - Full 3D→plane-stress material transformation
    - Diagonal symmetric BCs (ANSYS-compatible)
    - Multiple strain extraction methods
    """

    def __init__(self, params: Dict | IosipescuParams):
        """Initialize solver with parameters."""
        if not isinstance(params, IosipescuParams):
            params = IosipescuParams(**params)  # type: ignore
        self.p = params

        # Runtime containers (initialized during solve)
        self.nodes: np.ndarray = np.zeros((0, 2))
        self.elements: np.ndarray = np.zeros((0, 4), dtype=int)
        self.S2_elements: List[int] = []
        self.C_global: np.ndarray | None = None
        self.displacement: np.ndarray | None = None

    # =========================================================================
    # MATERIAL MODEL
    # =========================================================================
    def compute_stiffness_matrix_transformed(self) -> np.ndarray:
        """
        Compute in-plane stiffness matrix C (3×3) for plane stress.

        For plane stress in the LR material plane (T out-of-plane):
          1. Build 2D stiffness in material frame (L,R coordinates)
          2. Rotate to global frame using 2D transformation

        This matches ANSYS TB,ANEL + PLANE182 behavior exactly.

        Returns:
            C_plane: 3×3 stiffness [σ] = C [ε] with engineering shear γxy
                     Order: [xx, yy, xy]
        """
        p = self.p

        # --- Step 1: Build 2D stiffness in LR material plane ---
        # Symmetric Poisson's ratio
        nu_RL = (p.ER / p.EL) * p.CPLR

        # Common denominator for plane stress
        denom = 1.0 - p.CPLR * nu_RL

        # Reduced stiffness matrix in material frame (L, R, LR)
        # Qbar in material coordinates [L, R, LR]
        Q_mat = np.zeros((3, 3))
        Q_mat[0, 0] = p.EL / denom  # Q_LL
        Q_mat[1, 1] = p.ER / denom  # Q_RR
        Q_mat[0, 1] = p.CPLR * p.ER / denom  # Q_LR
        Q_mat[1, 0] = Q_mat[0, 1]  # Symmetry
        Q_mat[2, 2] = p.GLR  # Q_LR,LR (shear)

        # --- Step 2: Rotate to global frame (2D transformation) ---
        th = np.deg2rad(p.ang)
        c, s = np.cos(th), np.sin(th)
        c2, s2 = c ** 2, s ** 2
        cs = c * s

        # Transformation matrix for plane stress (Reuter matrix for engineering strain)
        # Transforms [εL, εR, γLR]^T to [εx, εy, γxy]^T
        T = np.array([
            [c2, s2, 2 * cs],
            [s2, c2, -2 * cs],
            [-cs, cs, c2 - s2]
        ])

        # Transformed stiffness: Q' = T^(-T) @ Q @ T^(-1)
        # For orthogonal transformations in 2D: T^(-1) = T^T (with proper scaling)
        # But for Reuter matrix, we use: Q' = inv(T)^T @ Q @ inv(T)
        T_inv = np.linalg.inv(T)
        self.C_global = T_inv.T @ Q_mat @ T_inv

        # Report
        Q = self.C_global / 1e3  # Convert to GPa for display
        print(f"Material stiffness (θ={p.ang}°):")
        print(f"  Q11={Q[0, 0]:.3f} GPa, Q22={Q[1, 1]:.3f} GPa")
        print(f"  Q12={Q[0, 1]:.3f} GPa, Q66={Q[2, 2]:.3f} GPa")

        return self.C_global

    # =========================================================================
    # MESH GENERATION
    # =========================================================================
    def generate_mesh(self) -> Tuple[int, int]:
        """
        Generate structured quad mesh for Iosipescu specimen.

        Returns:
            (n_nodes, n_elements): Mesh statistics
        """
        p = self.p

        # Define 5 x-regions with divisions
        x_regions = [
            (-p.IBC, -(p.IBC - p.OBC), p.LdH1),  # Left outer
            (-(p.IBC - p.OBC), 0.0, p.LdH2),  # Left transition
            (0.0, p.L, p.LdH3),  # S2 gauge section
            (p.L, p.L + (p.IBC - p.OBC), p.LdH2),  # Right transition
            (p.L + (p.IBC - p.OBC), p.L + p.IBC, p.LdH1)  # Right outer
        ]

        # Create x-coordinates (avoiding duplicates at boundaries)
        xs: List[float] = []
        for k, (xa, xb, n) in enumerate(x_regions):
            seg = np.linspace(xa, xb, n + 1)
            xs.extend(seg if k == 0 else seg[1:])
        xs_arr = np.array(xs)

        # Create y-coordinates
        ys = np.linspace(-p.W / 2.0, p.W / 2.0, p.LdV + 1)

        # Generate nodes in structured grid
        nodes = []
        idmap: Dict[Tuple[int, int], int] = {}
        nid = 0
        for j, y in enumerate(ys):
            for i, x in enumerate(xs_arr):
                nodes.append([x, y])
                idmap[(i, j)] = nid
                nid += 1
        self.nodes = np.array(nodes)

        # Generate quad elements (CCW node ordering)
        nx, ny = len(xs_arr), len(ys)
        elems: List[List[int]] = []
        for j in range(ny - 1):
            for i in range(nx - 1):
                n1 = idmap[(i, j)]
                n2 = idmap[(i + 1, j)]
                n3 = idmap[(i + 1, j + 1)]
                n4 = idmap[(i, j + 1)]
                elems.append([n1, n2, n3, n4])
        self.elements = np.array(elems, dtype=int)

        # Identify S2 elements (gauge section: x ∈ [0, L])
        x0_idx = p.LdH1 + p.LdH2
        x1_idx = x0_idx + p.LdH3
        self.S2_elements = []
        stride = nx - 1
        for j in range(ny - 1):
            for i in range(x0_idx, x1_idx):
                idx = j * stride + i
                if idx < len(self.elements):
                    self.S2_elements.append(idx)

        print(f"Mesh generated: {len(self.nodes)} nodes, {len(self.elements)} elements")
        print(f"  S2 region: {len(self.S2_elements)} elements")

        return len(self.nodes), len(self.elements)

    # =========================================================================
    # ELEMENT OPERATIONS
    # =========================================================================
    @staticmethod
    def _quad_B_matrix(elem_nodes: np.ndarray, xi: float, eta: float) -> Tuple[np.ndarray, float]:
        """
        Compute strain-displacement matrix B and Jacobian determinant.

        Args:
            elem_nodes: (4, 2) array of element node coordinates
            xi, eta: Isoparametric coordinates ∈ [-1, 1]

        Returns:
            B: (3, 8) strain-displacement matrix for [εx, εy, γxy]
            detJ: Jacobian determinant
        """
        # Shape function derivatives in isoparametric space
        dN_xi = 0.25 * np.array([
            [-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)],
            [-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)]
        ])

        # Jacobian matrix and determinant
        J = dN_xi @ elem_nodes
        detJ = np.linalg.det(J)

        if detJ <= JACOBIAN_MIN:
            return np.zeros((3, 8)), detJ

        # Shape function derivatives in physical space
        dN_xy = np.linalg.inv(J) @ dN_xi

        # Assemble B matrix: ε = B u, with u = [u1x, u1y, u2x, u2y, ...]
        B = np.zeros((3, 8))
        for i in range(4):
            B[0, 2 * i] = dN_xy[0, i]  # εx = ∂ux/∂x
            B[1, 2 * i + 1] = dN_xy[1, i]  # εy = ∂uy/∂y
            B[2, 2 * i] = dN_xy[1, i]  # γxy = ∂ux/∂y + ∂uy/∂x
            B[2, 2 * i + 1] = dN_xy[0, i]

        return B, detJ

    def element_stiffness(self, elem_nodes: np.ndarray) -> Tuple[np.ndarray, int]:
        """
        Compute element stiffness matrix using 2×2 Gauss integration.

        Args:
            elem_nodes: (4, 2) array of element node coordinates

        Returns:
            ke: (8, 8) element stiffness matrix
            bad_gp: Number of integration points with bad Jacobian
        """
        assert self.C_global is not None, "Call compute_stiffness_matrix_transformed() first"

        g = GAUSS_COORD
        gauss_points = [(-g, -g), (g, -g), (g, g), (-g, g)]

        ke = np.zeros((8, 8))
        bad_gp = 0

        for xi, eta in gauss_points:
            B, detJ = self._quad_B_matrix(elem_nodes, xi, eta)
            if detJ <= JACOBIAN_MIN:
                bad_gp += 1
                continue
            # ke = ∫ B^T C B dV = ∑_gp B^T C B detJ w t (w=1 for 2×2 Gauss)
            ke += B.T @ self.C_global @ B * detJ * self.p.Th

        return ke, bad_gp

    # =========================================================================
    # ASSEMBLY & BOUNDARY CONDITIONS
    # =========================================================================
    def assemble_global_stiffness(self) -> lil_matrix:
        """
        Assemble global stiffness matrix K.

        Returns:
            K: (2N, 2N) sparse stiffness matrix
        """
        n_dof = 2 * len(self.nodes)
        K = lil_matrix((n_dof, n_dof))

        total_bad, bad_elems = 0, 0
        for elem in self.elements:
            ke, bad = self.element_stiffness(self.nodes[elem])
            if bad:
                total_bad += bad
                bad_elems += 1

            # Assemble into global matrix
            dof = np.zeros(8, dtype=int)
            for i, n in enumerate(elem):
                dof[2 * i] = 2 * n
                dof[2 * i + 1] = 2 * n + 1

            for i in range(8):
                K[dof[i], dof] += ke[i, :]

        if total_bad:
            ratio = 100.0 * bad_elems / len(self.elements)
            print(f"  ⚠ Warning: {total_bad} bad Gauss points in {bad_elems} elements ({ratio:.1f}%)")

        return K

    def apply_boundary_conditions(self) -> Tuple[List[int], List[float]]:
        """
        Apply diagonal symmetric grip boundary conditions (ANSYS-compatible).

        Grip pattern:
          LEFT GRIPS (uy = 0):
            - Upper (SHORT): y = +W/2, x ∈ [-IBC, -(IBC-OBC)]
            - Lower (LONG):  y = -W/2, x ∈ [-IBC, 0]

          RIGHT GRIPS (uy = -uy_prescribed):
            - Upper (LONG):  y = +W/2, x ∈ [L, L+IBC]
            - Lower (SHORT): y = -W/2, x ∈ [L+(IBC-OBC), L+IBC]

        Additionally: Fix one left-upper node in x to prevent rigid motion.

        Returns:
            bc_dof: List of constrained DOF indices
            bc_val: List of prescribed values
        """
        p = self.p
        tol = TOLERANCE

        bc_dof: List[int] = []
        bc_val: List[float] = []

        # Track nodes for reporting
        left_upper, left_lower = [], []
        right_upper, right_lower = [], []

        for i, (x, y) in enumerate(self.nodes):
            # Left upper grip (SHORT, uy=0)
            if abs(y - p.W / 2) < tol and (-p.IBC - tol) <= x <= (-(p.IBC - p.OBC) + tol):
                bc_dof.append(2 * i + 1)
                bc_val.append(0.0)
                left_upper.append(i)

            # Left lower grip (LONG, uy=0)
            elif abs(y + p.W / 2) < tol and (-p.IBC - tol) <= x <= (0.0 + tol):
                bc_dof.append(2 * i + 1)
                bc_val.append(0.0)
                left_lower.append(i)

            # Right upper grip (LONG, uy=-uy)
            elif abs(y - p.W / 2) < tol and (p.L - tol) <= x <= (p.L + p.IBC + tol):
                bc_dof.append(2 * i + 1)
                bc_val.append(-p.uy)
                right_upper.append(i)

            # Right lower grip (SHORT, uy=-uy)
            elif abs(y + p.W / 2) < tol and (p.L + (p.IBC - p.OBC) - tol) <= x <= (p.L + p.IBC + tol):
                bc_dof.append(2 * i + 1)
                bc_val.append(-p.uy)
                right_lower.append(i)

        # Fix one left-upper node in x (prevent rigid translation)
        if left_upper:
            leftmost = min(left_upper, key=lambda i: self.nodes[i][0])
            bc_dof.append(2 * leftmost)
            bc_val.append(0.0)

        # Report
        print("Boundary conditions applied:")
        print(f"  Left upper (uy=0):      {len(left_upper)} nodes")
        print(f"  Left lower (uy=0):      {len(left_lower)} nodes")
        print(f"  Right upper (uy=-{p.uy}): {len(right_upper)} nodes")
        print(f"  Right lower (uy=-{p.uy}): {len(right_lower)} nodes")

        if len(set(d // 2 for d in bc_dof)) < 3:
            warnings.warn("Possibly insufficient constraints - check for rigid modes")

        return bc_dof, bc_val

    # =========================================================================
    # SOLVER
    # =========================================================================
    def solve(self) -> np.ndarray:
        """
        Solve the FEM problem: K U = F.

        Uses exact Dirichlet enforcement (DOF elimination, not penalty).

        Returns:
            U: (2N,) displacement vector
        """
        print("=" * 60)
        print("STARTING FEM ANALYSIS")
        print("=" * 60)

        # Step 1: Material
        print("\n[1/5] Computing material stiffness...")
        self.compute_stiffness_matrix_transformed()

        # Step 2: Mesh
        print("\n[2/5] Generating mesh...")
        self.generate_mesh()

        # Step 3: Assembly
        print("\n[3/5] Assembling global stiffness...")
        K = self.assemble_global_stiffness().tocsr()
        print(f"  K: {K.shape[0]} DOFs, {K.nnz:,} non-zeros")

        # Step 4: Boundary conditions
        print("\n[4/5] Applying boundary conditions...")
        bc_dof, bc_val = self.apply_boundary_conditions()

        # Step 5: Solve
        print("\n[5/5] Solving linear system...")
        n_dof = 2 * len(self.nodes)
        U = np.zeros(n_dof)

        # Mark constrained DOFs
        constrained = np.zeros(n_dof, dtype=bool)
        vals = np.zeros(n_dof)
        for d, v in zip(bc_dof, bc_val):
            constrained[d] = True
            vals[d] = v
        free = ~constrained

        # External forces (zero for this problem)
        F = np.zeros(n_dof)

        # Partitioned system: K_ff U_f = F_f - K_fc U_c
        Kff = K[free][:, free]
        Kfc = K[free][:, constrained]
        Ff = F[free]

        rhs = Ff - Kfc @ vals[constrained]

        try:
            Uf = spsolve(Kff, rhs)
        except RuntimeError as e:
            raise RuntimeError(
                f"Solver failed! Check:\n"
                f"  - Material properties (all positive?)\n"
                f"  - Boundary conditions (sufficient constraints?)\n"
                f"  - Mesh quality (badly distorted elements?)\n"
                f"Original error: {e}"
            )

        # Assemble full solution
        U[constrained] = vals[constrained]
        U[free] = Uf

        if not np.isfinite(U).all():
            raise ValueError("Solution contains non-finite values")

        self.displacement = U

        print(f"  ✓ Solution complete")
        print(f"  Max displacement: |U| = {np.abs(U).max():.3e} mm")
        print("=" * 60)

        return U

    # =========================================================================
    # POST-PROCESSING: STRAIN EXTRACTION
    # =========================================================================
    @staticmethod
    def _quad_area(elem_nodes: np.ndarray) -> float:
        """Compute quadrilateral area using shoelace formula."""
        x, y = elem_nodes[:, 0], elem_nodes[:, 1]
        return 0.5 * abs(
            x[0] * y[1] - x[1] * y[0] +
            x[1] * y[2] - x[2] * y[1] +
            x[2] * y[3] - x[3] * y[2] +
            x[3] * y[0] - x[0] * y[3]
        )

    def _elem_displacement_vector(self, elem: np.ndarray) -> np.ndarray:
        """Extract element displacement vector [u1x, u1y, u2x, u2y, ...]."""
        assert self.displacement is not None
        ue = np.zeros(8)
        for i, n in enumerate(elem):
            ue[2 * i] = self.displacement[2 * n]
            ue[2 * i + 1] = self.displacement[2 * n + 1]
        return ue

    def compute_nodal_averaged_strains(self) -> np.ndarray:
        """
        Compute nodal-averaged engineering strains using integration-weighted method.

        This is the ANSYS-like approach:
        - Evaluate strains at Gauss points
        - Distribute to nodes using shape function weights
        - Average contributions from all connected elements

        More stable than extrapolation (no overshoot).

        Returns:
            eps_node: (n_nodes, 3) array of [εx, εy, γxy] at each node
        """
        assert self.C_global is not None and self.displacement is not None

        nnode = len(self.nodes)
        numerator = np.zeros((nnode, 3))
        denominator = np.zeros(nnode)

        g = GAUSS_COORD
        gauss_points = [(-g, -g), (g, -g), (g, g), (-g, g)]
        weights = [1.0, 1.0, 1.0, 1.0]  # Weight for each GP

        for elem in self.elements:
            xy = self.nodes[elem]
            ue = self._elem_displacement_vector(elem)

            for (xi, eta), w in zip(gauss_points, weights):
                # Shape functions at this Gauss point
                N = 0.25 * np.array([
                    (1 - xi) * (1 - eta),
                    (1 + xi) * (1 - eta),
                    (1 + xi) * (1 + eta),
                    (1 - xi) * (1 + eta)
                ])

                # Strain at GP
                B, detJ = self._quad_B_matrix(xy, xi, eta)
                if detJ <= JACOBIAN_MIN:
                    continue

                eps = B @ ue  # Engineering strains [εx, εy, γxy]

                # Integration weight (area measure)
                wA = detJ * w

                # Distribute to connected nodes
                for k, n in enumerate(elem):
                    numerator[n] += N[k] * wA * eps
                    denominator[n] += N[k] * wA

        # Avoid division by zero
        denominator[denominator == 0] = 1.0
        eps_node = numerator / denominator[:, None]

        return eps_node

    def compute_s2_results(
            self,
            method: Literal["element_average", "centroid", "nodal_average", "gauss_points"] = "element_average"
    ) -> List[Dict[str, float]]:
        """
        Extract S2 region results using specified method.

        Methods:
          - "element_average": Average of 4 Gauss points (ANSYS ETABLE, EPEL)
          - "centroid": Single evaluation at element center (ξ=η=0)
          - "nodal_average": Average of 4 corner nodal values
          - "gauss_points": One row per Gauss point (4 rows/element)

        Returns:
            List of dictionaries with keys:
              area, xc, yc, ux, uy, ex, ey, exy, [sx, sy, sxy if computed]
        """
        assert self.C_global is not None and self.displacement is not None

        if method == "element_average":
            return self._compute_element_average()
        elif method == "centroid":
            return self._compute_centroid()
        elif method == "nodal_average":
            return self._compute_nodal_average()
        elif method == "gauss_points":
            return self._compute_gauss_points()
        else:
            raise ValueError(f"Unknown method: {method}")

    def _compute_element_average(self) -> List[Dict[str, float]]:
        """Average strains/stresses over 4 Gauss points (ANSYS EPEL-like)."""
        results = []
        g = GAUSS_COORD
        gauss_points = [(-g, -g), (g, -g), (g, g), (-g, g)]

        for eidx in self.S2_elements:
            elem = self.elements[eidx]
            xy = self.nodes[elem]
            A = self._quad_area(xy)
            if A < JACOBIAN_MIN:
                continue

            ue = self._elem_displacement_vector(elem)

            # Average displacement (nodal average)
            ux = float((ue[0] + ue[2] + ue[4] + ue[6]) / 4.0)
            uy = float((ue[1] + ue[3] + ue[5] + ue[7]) / 4.0)

            # Average strains over GPs
            eps_sum = np.zeros(3)
            sig_sum = np.zeros(3)
            valid = 0

            for xi, eta in gauss_points:
                B, detJ = self._quad_B_matrix(xy, xi, eta)
                if detJ <= JACOBIAN_MIN:
                    continue

                eps = B @ ue
                sig = self.C_global @ eps
                eps_sum += eps
                sig_sum += sig
                valid += 1

            if valid == 0:
                continue

            eps_avg = eps_sum / valid
            sig_avg = sig_sum / valid

            xc = float(xy[:, 0].mean())
            yc = float(xy[:, 1].mean())

            results.append({
                "area": A, "xc": xc, "yc": yc,
                "ux": ux, "uy": uy,
                "ex": float(eps_avg[0]), "ey": float(eps_avg[1]), "exy": float(eps_avg[2]),
                "sx": float(sig_avg[0]), "sy": float(sig_avg[1]), "sxy": float(sig_avg[2]),
            })

        return results

    def _compute_centroid(self) -> List[Dict[str, float]]:
        """Evaluate at element centroid (ξ=η=0)."""
        results = []

        for eidx in self.S2_elements:
            elem = self.elements[eidx]
            xy = self.nodes[elem]
            A = self._quad_area(xy)
            if A < JACOBIAN_MIN:
                continue

            ue = self._elem_displacement_vector(elem)

            # Centroid evaluation
            B, detJ = self._quad_B_matrix(xy, 0.0, 0.0)
            if detJ <= JACOBIAN_MIN:
                continue

            eps = B @ ue
            sig = self.C_global @ eps

            xc = float(xy[:, 0].mean())
            yc = float(xy[:, 1].mean())
            ux = float((ue[0] + ue[2] + ue[4] + ue[6]) / 4.0)
            uy = float((ue[1] + ue[3] + ue[5] + ue[7]) / 4.0)

            results.append({
                "area": A, "xc": xc, "yc": yc,
                "ux": ux, "uy": uy,
                "ex": float(eps[0]), "ey": float(eps[1]), "exy": float(eps[2]),
                "sx": float(sig[0]), "sy": float(sig[1]), "sxy": float(sig[2]),
            })

        return results

    def _compute_nodal_average(self) -> List[Dict[str, float]]:
        """Average of 4 corner nodal strains."""
        eps_node = self.compute_nodal_averaged_strains()
        results = []

        for eidx in self.S2_elements:
            elem = self.elements[eidx]
            xy = self.nodes[elem]
            A = self._quad_area(xy)
            if A < JACOBIAN_MIN:
                continue

            ue = self._elem_displacement_vector(elem)

            # Average nodal strains
            eps_el = eps_node[elem].mean(axis=0)
            sig_el = self.C_global @ eps_el

            xc = float(xy[:, 0].mean())
            yc = float(xy[:, 1].mean())
            ux = float((ue[0] + ue[2] + ue[4] + ue[6]) / 4.0)
            uy = float((ue[1] + ue[3] + ue[5] + ue[7]) / 4.0)

            results.append({
                "area": A, "xc": xc, "yc": yc,
                "ux": ux, "uy": uy,
                "ex": float(eps_el[0]), "ey": float(eps_el[1]), "exy": float(eps_el[2]),
                "sx": float(sig_el[0]), "sy": float(sig_el[1]), "sxy": float(sig_el[2]),
            })

        return results

    def _compute_gauss_points(self) -> List[Dict[str, float]]:
        """One row per Gauss point (4 rows/element, Area=A/4)."""
        results = []
        g = GAUSS_COORD
        gauss_points = [(-g, -g), (g, -g), (g, g), (-g, g)]

        for eidx in self.S2_elements:
            elem = self.elements[eidx]
            xy = self.nodes[elem]
            A = self._quad_area(xy)
            if A < JACOBIAN_MIN:
                continue

            ue = self._elem_displacement_vector(elem)

            for xi, eta in gauss_points:
                # Shape functions for interpolation
                N = 0.25 * np.array([
                    (1 - xi) * (1 - eta),
                    (1 + xi) * (1 - eta),
                    (1 + xi) * (1 + eta),
                    (1 - xi) * (1 + eta)
                ])

                xg = float(N @ xy[:, 0])
                yg = float(N @ xy[:, 1])
                ux = float(sum(N[i] * ue[2 * i] for i in range(4)))
                uy = float(sum(N[i] * ue[2 * i + 1] for i in range(4)))

                B, detJ = self._quad_B_matrix(xy, xi, eta)
                if detJ <= JACOBIAN_MIN:
                    continue

                eps = B @ ue
                sig = self.C_global @ eps

                results.append({
                    "area": A / 4.0, "xc": xg, "yc": yg,
                    "ux": ux, "uy": uy,
                    "ex": float(eps[0]), "ey": float(eps[1]), "exy": float(eps[2]),
                    "sx": float(sig[0]), "sy": float(sig[1]), "sxy": float(sig[2]),
                })

        return results

    def compute_reaction_force(self) -> float:
        """
        Compute total reaction force on right grips.

        Returns:
            P: Total force magnitude (N)
        """
        assert self.displacement is not None

        K = self.assemble_global_stiffness().tocsr()
        F_react = K @ self.displacement

        p = self.p
        tol = TOLERANCE
        total = 0.0

        for i, (x, y) in enumerate(self.nodes):
            # Right upper grip
            if (p.L - tol) <= x <= (p.L + p.IBC + tol) and abs(y - p.W / 2) < tol:
                total += F_react[2 * i + 1]
            # Right lower grip
            if (p.L + (p.IBC - p.OBC) - tol) <= x <= (p.L + p.IBC + tol) and abs(y + p.W / 2) < tol:
                total += F_react[2 * i + 1]

        return abs(float(total))

    # =========================================================================
    # EXPORT
    # =========================================================================
    def export_results(
            self,
            method: Literal["element_average", "centroid", "nodal_average", "gauss_points"] = "element_average",
            scalar_file: str = "scalarsFE.csv",
            field_file: str = "FEM2VFM.csv"
    ) -> None:
        """
        Export results to CSV files for VFM processing.

        Args:
            method: Strain extraction method
            scalar_file: Filename for global results
            field_file: Filename for field data
        """
        # Compute reaction force
        P = self.compute_reaction_force()

        # Export scalars
        print(f"\n{'=' * 60}")
        print("RESULTS SUMMARY")
        print(f"{'=' * 60}")
        print(f"Angle  = {self.p.ang:.3f}°")
        print(f"Length = {self.p.L:.3f} mm")
        print(f"Width  = {self.p.W:.3f} mm")
        print(f"Thick  = {self.p.Th:.3f} mm")
        print(f"Force  = {-P:.3f} N")
        print(f"{'=' * 60}")

        with open(scalar_file, "w") as f:
            f.write(f"Angle = {self.p.ang:.3f}\n")
            f.write(f"Length = {self.p.L:.3f}\n")
            f.write(f"Width = {self.p.W:.3f}\n")
            f.write(f"Thick = {self.p.Th:.3f}\n")
            f.write(f"P = {-P:.3f}\n")

        # Export field data
        rows = self.compute_s2_results(method=method)

        with open(field_file, "w") as f:
            f.write("Area,X_Coord,Y_Coord,U_X,U_Y,Eps_X,Eps_Y,Eps_XY\n")
            for r in rows:
                f.write(
                    f"{r['area']:.4f},{r['xc']:.4f},{r['yc']:.4f},"
                    f"{r['ux']:.9E},{r['uy']:.9E},"
                    f"{r['ex']:.9E},{r['ey']:.9E},{r['exy']:.9E}\n"
                )

        print(f"\n✓ Exported {scalar_file} (global results)")
        print(f"✓ Exported {field_file} ({len(rows)} rows, method='{method}')")

    # =========================================================================
    # VISUALIZATION
    # =========================================================================
    def plot_mesh(self, show_bc: bool = True, figsize: Tuple[int, int] = (12, 6)):
        """Plot mesh with optional boundary conditions."""
        fig, ax = plt.subplots(figsize=figsize)

        # Elements
        for elem in self.elements:
            coords = np.vstack([self.nodes[elem], self.nodes[elem][0]])
            ax.plot(coords[:, 0], coords[:, 1], "k-", lw=0.5, alpha=0.6)

        # Nodes
        ax.plot(self.nodes[:, 0], self.nodes[:, 1], "b.", ms=2, alpha=0.5)

        # S2 highlight
        for eidx in self.S2_elements:
            poly = Polygon(self.nodes[self.elements[eidx]],
                           facecolor="lightblue", edgecolor="blue",
                           alpha=0.25, lw=0.5)
            ax.add_patch(poly)

        if show_bc:
            p = self.p
            tol = TOLERANCE

            # Helper to filter nodes
            def filter_nodes(condition):
                return np.array([n for n in self.nodes if condition(*n)])

            # BC regions
            lu = filter_nodes(lambda x, y: abs(y - p.W / 2) < tol and -p.IBC <= x <= -(p.IBC - p.OBC))
            ll = filter_nodes(lambda x, y: abs(y + p.W / 2) < tol and -p.IBC <= x <= 0)
            ru = filter_nodes(lambda x, y: abs(y - p.W / 2) < tol and p.L <= x <= p.L + p.IBC)
            rl = filter_nodes(lambda x, y: abs(y + p.W / 2) < tol and p.L + (p.IBC - p.OBC) <= x <= p.L + p.IBC)

            # Plot BC nodes
            handles = []
            for pts, style, label in [
                (lu, "go", "Left Upper (uy=0)"),
                (ll, "cs", "Left Lower (uy=0)"),
                (ru, "r^", f"Right Upper (uy=-{p.uy})"),
                (rl, "mv", f"Right Lower (uy=-{p.uy})"),
            ]:
                if len(pts):
                    h, = ax.plot(pts[:, 0], pts[:, 1], style, ms=6, zorder=10, label=label)
                    handles.append(h)

            # Reference lines
            for x in [-p.IBC, -(p.IBC - p.OBC), p.L, p.L + (p.IBC - p.OBC), p.L + p.IBC]:
                ax.axvline(x, ls=":", alpha=0.5, color="gray")

        # S2 box
        ax.plot([0, p.L, p.L, 0, 0],
                [-p.W / 2, -p.W / 2, p.W / 2, p.W / 2, -p.W / 2],
                "b--", lw=1.5, alpha=0.6, label="S2 region")

        ax.set_xlabel("X (mm)")
        ax.set_ylabel("Y (mm)")
        ax.set_title(f"Iosipescu Mesh (θ={self.p.ang}°, OBC={self.p.OBC}mm)")
        ax.axis("equal")
        ax.grid(True, alpha=0.3)
        if show_bc:
            ax.legend(fontsize=8, loc="best")

        plt.tight_layout()
        return fig

    def plot_deformed(self, scale: float = 10.0, figsize: Tuple[int, int] = (12, 6)):
        """Plot undeformed and deformed mesh."""
        assert self.displacement is not None

        fig, ax = plt.subplots(figsize=figsize)

        # Deformed coordinates
        nodes_def = self.nodes.copy()
        for i in range(len(self.nodes)):
            nodes_def[i, 0] += self.displacement[2 * i] * scale
            nodes_def[i, 1] += self.displacement[2 * i + 1] * scale

        # Plot
        first = True
        for elem in self.elements:
            c0 = np.vstack([self.nodes[elem], self.nodes[elem][0]])
            c1 = np.vstack([nodes_def[elem], nodes_def[elem][0]])

            ax.plot(c0[:, 0], c0[:, 1], "b-", lw=0.5, alpha=0.4,
                    label="Undeformed" if first else "")
            ax.plot(c1[:, 0], c1[:, 1], "r-", lw=0.8, alpha=0.7,
                    label="Deformed" if first else "")
            first = False

        ax.set_xlabel("X (mm)")
        ax.set_ylabel("Y (mm)")
        ax.set_title(f"Deformed Mesh (scale={scale}×)")
        ax.axis("equal")
        ax.grid(True, alpha=0.3)
        ax.legend()
        plt.tight_layout()
        return fig

    def plot_field(self, field: str, figsize: Tuple[int, int] = (10, 8), cmap: str = "jet"):
        """
        Plot contour of field quantity over S2 region.

        Args:
            field: One of {ux, uy, ex, ey, exy, sx, sy, sxy}
            figsize: Figure size
            cmap: Colormap name
        """
        # Get data at centroids
        data = self._compute_centroid()
        if not data:
            warnings.warn(f"No valid data for field '{field}'")
            return None

        x = [d["xc"] for d in data]
        y = [d["yc"] for d in data]
        z = [d[field] for d in data]

        fig, ax = plt.subplots(figsize=figsize)

        # Triangulated contour
        tri = Triangulation(x, y)
        cf = ax.tricontourf(tri, z, levels=20, cmap=cmap)
        ax.tricontour(tri, z, levels=20, colors="k", lw=0.4, alpha=0.3)

        cbar = fig.colorbar(cf, ax=ax)

        # Labels
        labels = {
            "ux": ("Displacement u_x", "mm"),
            "uy": ("Displacement u_y", "mm"),
            "ex": ("Strain ε_x", ""),
            "ey": ("Strain ε_y", ""),
            "exy": ("Shear strain γ_xy", ""),
            "sx": ("Stress σ_x", "MPa"),
            "sy": ("Stress σ_y", "MPa"),
            "sxy": ("Shear stress τ_xy", "MPa"),
        }
        label, unit = labels.get(field, (field, ""))
        cbar.set_label(f"{label} ({unit})" if unit else label)

        # S2 boundary
        p = self.p
        ax.plot([0, p.L, p.L, 0, 0],
                [-p.W / 2, -p.W / 2, p.W / 2, p.W / 2, -p.W / 2],
                "r--", lw=1.2, alpha=0.6)

        ax.set_xlabel("X (mm)")
        ax.set_ylabel("Y (mm)")
        ax.set_title(f"{label} (S2 region)")
        ax.axis("equal")
        ax.grid(True, alpha=0.25)
        plt.tight_layout()
        return fig

    def save_all_plots(self, dpi: int = 150, prefix: str = "") -> None:
        """Save all standard plots to files."""
        print("\nGenerating plots...")

        plots = [
            ("mesh", lambda: self.plot_mesh()),
            ("deformed", lambda: self.plot_deformed(scale=10.0)),
        ]

        fields = [
            ("ux", "RdBu_r"), ("uy", "RdBu_r"),
            ("ex", "coolwarm"), ("ey", "coolwarm"), ("exy", "coolwarm"),
            ("sx", "seismic"), ("sy", "seismic"), ("sxy", "seismic"),
        ]

        for i, (name, func) in enumerate(plots, 1):
            fig = func()
            filename = f"{prefix}{i:02d}_{name}.png"
            fig.savefig(filename, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
            print(f"  ✓ {filename}")

        for i, (field, cmap) in enumerate(fields, len(plots) + 1):
            fig = self.plot_field(field, cmap=cmap)
            if fig:
                filename = f"{prefix}{i:02d}_{field}_S2.png"
                fig.savefig(filename, dpi=dpi, bbox_inches="tight")
                plt.close(fig)
                print(f"  ✓ {filename}")

        print("All plots saved.")


# =============================================================================
# MAIN
# =============================================================================
def main():
    """Example usage."""

    # Define test parameters
    params = IosipescuParams(
        # Geometry (mm)
        ang=0.0,
        L=34.0,
        W=20.0,
        Th=5.0,
        tL=80.0,
        OBC=15.0,
        IBC=23.0,

        # Material (MPa)
        EL=15100.0,
        ER=1910.0,
        ET=1010.0,
        CPLR=0.470,
        CPTL=0.051,
        CPRT=0.586,
        GLR=1109.0,
        GLT=1096.0,
        GRT=176.0,

        # Mesh
        LdH3=85,
        LdV=50,
        LdH1=38,
        LdH2=20,

        # Loading (mm)
        uy=1.0,
    )

    # Create and run solver
    fem = IosipescuFEM(params)
    fem.solve()

    # Export results (try different methods)
    print("\n" + "=" * 60)
    print("EXPORTING RESULTS")
    print("=" * 60)

    fem.export_results(method="element_average")

    # Optionally try other methods:
    # fem.export_results(method="centroid", field_file="FEM2VFM_centroid.csv")
    # fem.export_results(method="nodal_average", field_file="FEM2VFM_nodal.csv")

    # Generate plots
    fem.save_all_plots(dpi=150)

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print("Files generated:")
    print("  • scalarsFE.csv (global results)")
    print("  • FEM2VFM.csv (field data for VFM)")
    print("  • 01_mesh.png ... 10_sxy_S2.png (plots)")
    print("=" * 60)


if __name__ == "__main__":
    main()