"""
Virtual Fields Method: Three sets of virtual fields with detailed validation

Based on:
"The Virtual Fields Methods: Extracting constitutive mechanical parameters
from full-field deformation measurements" by F. Pierron, M. Grédiac
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import io


class IosipescuVFM:
    """
    Debug version of VFM analysis with three virtual field sets
    """

    def __init__(self, data_source='csv'):
        self.data_source = data_source
        self.data = {}
        self.results = {}

    def load_data_from_csv(self, csv_dir='.'):
        # Load scalar parameters from scalarsFE.csv
        scalar_path = os.path.join(csv_dir, 'scalarsFE.csv')

        if os.path.exists(scalar_path):
            with open(scalar_path, 'r') as f:
                for line in f:
                    if '=' in line:
                        key, value = line.strip().split('=')
                        key = key.strip()
                        value = float(value.strip())

                        # Map the scalar parameters
                        if key == 'Length':
                            self.data['L'] = value
                        elif key == 'Width':
                            self.data['w'] = value
                        elif key == 'Thick':
                            self.data['t'] = value
                        elif key == 'P':
                            self.data['F'] = value

        # Load FEM data from FEM2VFM.csv
        fem_path = os.path.join(csv_dir, 'FEM2VFM.csv')

        if os.path.exists(fem_path):
            # Read the space-delimited file
            fem_data = pd.read_csv(fem_path, sep=r'\s+')

            # Map the field data
            self.data['X1'] = fem_data['X_Coord'].values
            self.data['X2'] = fem_data['Y_Coord'].values
            self.data['Eps1'] = fem_data['Eps_X'].values
            self.data['Eps2'] = fem_data['Eps_Y'].values
            self.data['Eps6'] = fem_data['Eps_XY'].values

        # Validate data
        self._validate_data()

    def _validate_data(self):
        """Validate loaded data"""
        print("Data validation:")
        for key, value in self.data.items():
            if np.isscalar(value):
                print(f"  {key}: {value}")
            else:
                arr = np.array(value)
                print(f"  {key}: shape {arr.shape}, range [{arr.min():.3f}, {arr.max():.3f}]")

        # Check coordinate ranges
        X1, X2 = self.data['X1'], self.data['X2']
        L, w = self.data['L'], self.data['w']
        print(f"\nCoordinate validation:")
        print(f"  X1 range: [{X1.min():.3f}, {X1.max():.3f}], Expected: [0, {L}]")
        print(f"  X2 range: [{X2.min():.3f}, {X2.max():.3f}], Expected: [approx ±{w / 2}]")

    def load_data(self, source_path='.'):
        if self.data_source == 'csv':
            self.load_data_from_csv(source_path)
        else:
            raise ValueError("Only CSV loading implemented in debug version")

    def check_virtual_field_consistency(self):
        """Check if virtual fields satisfy boundary conditions"""
        X1 = self.data['X1']
        X2 = self.data['X2']
        L = self.data['L']

        print("\nVirtual field boundary condition checks:")

        # Virtual field 1: U1*=0, U2*=-x1
        # Should have U1*=U2*=0 at x1=0 and U1*=0, U2*=constant at x1=L
        u1_vf1 = np.zeros_like(X1)
        u2_vf1 = -X1

        # Check at x1=0
        mask_x1_0 = np.abs(X1) < 1e-6
        if np.any(mask_x1_0):
            print(f"  VF1 at x1=0: U1*={u1_vf1[mask_x1_0].mean():.6f}, U2*={u2_vf1[mask_x1_0].mean():.6f}")

        # Check at x1=L
        mask_x1_L = np.abs(X1 - L) < 1e-6
        if np.any(mask_x1_L):
            print(f"  VF1 at x1=L: U1*={u1_vf1[mask_x1_L].mean():.6f}, U2*={u2_vf1[mask_x1_L].mean():.6f}")

    def virtual_fields_set_1_debug(self):
        """First set of virtual fields with detailed debugging"""
        X1 = self.data['X1']
        X2 = self.data['X2']
        L = self.data['L']
        w = self.data['w']
        t = self.data['t']
        F = self.data['F']
        Eps1 = self.data['Eps1']
        Eps2 = self.data['Eps2']
        Eps6 = self.data['Eps6']

        print(f"\n" + "=" * 50)
        print("VIRTUAL FIELDS SET 1")
        print("=" * 50)
        print(f"Input parameters: L={L}, w={w}, t={t}, F={F}")

        A = np.zeros((4, 4))
        B = np.zeros(4)

        # VF1: U1*=0, U2*=-x1 → ε1*=0, ε2*=0, ε6*=-1
        A[0, 0] = 0  # Q11 coefficient
        A[0, 1] = 0  # Q22 coefficient
        A[0, 2] = 0  # Q12 coefficient
        A[0, 3] = np.mean(Eps6)  # Q66 coefficient
        B[0] = F / (w * t)
        print(np.mean(Eps6), F, w, t)

        print(f"VF1: U1*=0, U2*=-x1")
        print(f"     A[0,:] = {A[0, :]}")
        print(f"     B[0] = {B[0]}")

        # VF2: U1*=x1(L-x1)x2, U2*=x1^3/3-Lx1^2/2 → ε1*=(L-2x1)x2, ε2*=0, ε6*=0
        vf2_eps1_star = (L - 2 * X1) * X2
        A[1, 0] = -np.mean(Eps1 * vf2_eps1_star)
        A[1, 1] = 0
        A[1, 2] = -np.mean(Eps2 * vf2_eps1_star)
        A[1, 3] = 0
        B[1] = F * L * L / (6 * w * t)

        print(f"VF2: U1*=x1(L-x1)x2, U2*=x1³/3-Lx1²/2")
        print(f"     ε1* range = [{vf2_eps1_star.min():.3f}, {vf2_eps1_star.max():.3f}]")
        print(f"     A[1,:] = {A[1, :]}")
        print(f"     B[1] = {B[1]}")

        # VF3: U1*=0, U2*=x1(L-x1)x2 → ε1*=0, ε2*=x1(L-x1), ε6*=(L-2x1)x2
        vf3_eps2_star = X1 * (L - X1)
        vf3_eps6_star = (L - 2 * X1) * X2
        A[2, 0] = 0
        A[2, 1] = -np.mean(Eps2 * vf3_eps2_star)
        A[2, 2] = -np.mean(Eps1 * vf3_eps2_star)
        A[2, 3] = -np.mean(Eps6 * vf3_eps6_star)
        B[2] = 0

        print(f"VF3: U1*=0, U2*=x1(L-x1)x2")
        print(f"     ε2* range = [{vf3_eps2_star.min():.3f}, {vf3_eps2_star.max():.3f}]")
        print(f"     ε6* range = [{vf3_eps6_star.min():.3f}, {vf3_eps6_star.max():.3f}]")
        print(f"     A[2,:] = {A[2, :]}")
        print(f"     B[2] = {B[2]}")

        # VF4: U1*=L*sin(2πx1/L)/(2π), U2*=0 → ε1*=cos(2πx1/L), ε2*=0, ε6*=0
        vf4_eps1_star = np.cos(2 * np.pi * X1 / L)
        A[3, 0] = np.mean(Eps1 * vf4_eps1_star)
        A[3, 1] = 0
        A[3, 2] = np.mean(Eps2 * vf4_eps1_star)
        A[3, 3] = 0
        B[3] = 0

        print(f"VF4: U1*=L*sin(2πx1/L)/(2π), U2*=0")
        print(f"     ε1* range = [{vf4_eps1_star.min():.3f}, {vf4_eps1_star.max():.3f}]")
        print(f"     A[3,:] = {A[3, :]}")
        print(f"     B[3] = {B[3]}")

        return self._solve_system(A, B, "Set 1")

    def virtual_fields_set_2_debug(self):
        """Second set of virtual fields with detailed debugging"""
        X1 = self.data['X1']
        X2 = self.data['X2']
        L = self.data['L']
        w = self.data['w']
        t = self.data['t']
        F = self.data['F']
        Eps1 = self.data['Eps1']
        Eps2 = self.data['Eps2']
        Eps6 = self.data['Eps6']

        print(f"\n" + "=" * 50)
        print("VIRTUAL FIELDS SET 2")
        print("=" * 50)

        A = np.zeros((4, 4))
        B = np.zeros(4)

        # VF1: U1*=0, U2*=-x1 → ε1*=0, ε2*=0, ε6*=-1
        A[0, 0] = 0
        A[0, 1] = 0
        A[0, 2] = 0
        A[0, 3] = np.mean(Eps6)
        B[0] = F / (w * t)

        print(f"VF1: U1*=0, U2*=-x1")
        print(f"     A[0,:] = {A[0, :]}")
        print(f"     B[0] = {B[0]}")

        # VF2: U1*=x1(L-x1)x2, U2*=x1^3/3-Lx1^2/2 → ε1*=(L-2x1)x2, ε2*=0, ε6*=0
        vf2_eps1_star = (L - 2 * X1) * X2
        A[1, 0] = -np.mean(Eps1 * vf2_eps1_star)
        A[1, 1] = 0
        A[1, 2] = -np.mean(Eps2 * vf2_eps1_star)
        A[1, 3] = 0
        B[1] = F * L * L / (6 * w * t)

        print(f"VF2: U1*=x1(L-x1)x2, U2*=x1³/3-Lx1²/2")
        print(f"     ε1* range = [{vf2_eps1_star.min():.3f}, {vf2_eps1_star.max():.3f}]")
        print(f"     A[1,:] = {A[1, :]}")
        print(f"     B[1] = {B[1]}")

        # VF3: U1*=0, U2*=x1(L-x1)x2 → ε1*=0, ε2*=x1(L-x1), ε6*=(L-2x1)x2
        vf3_eps2_star = X1 * (L - X1)
        vf3_eps6_star = (L - 2 * X1) * X2
        A[2, 0] = 0
        A[2, 1] = -np.mean(Eps2 * vf3_eps2_star)
        A[2, 2] = -np.mean(Eps1 * vf3_eps2_star)
        A[2, 3] = -np.mean(Eps6 * vf3_eps6_star)
        B[2] = 0

        print(f"VF3: U1*=0, U2*=x1(L-x1)x2")
        print(f"     ε2* range = [{vf3_eps2_star.min():.3f}, {vf3_eps2_star.max():.3f}]")
        print(f"     ε6* range = [{vf3_eps6_star.min():.3f}, {vf3_eps6_star.max():.3f}]")
        print(f"     A[2,:] = {A[2, :]}")
        print(f"     B[2] = {B[2]}")

        # VF4: U1*=0, U2*=x1(L-x1)x2^3 → ε1*=0, ε2*=3x1(L-x1)x2^2, ε6*=(L-2x1)x2^3
        vf4_eps2_star = 3 * X1 * (L - X1) * X2 * X2
        vf4_eps6_star = (L - 2 * X1) * X2 * X2 * X2
        A[3, 0] = 0
        A[3, 1] = -np.mean(Eps2 * vf4_eps2_star)
        A[3, 2] = -np.mean(Eps1 * vf4_eps2_star)
        A[3, 3] = -np.mean(Eps6 * vf4_eps6_star)
        B[3] = 0

        print(f"VF4: U1*=0, U2*=x1(L-x1)x2³")
        print(f"     ε2* range = [{vf4_eps2_star.min():.3f}, {vf4_eps2_star.max():.3f}]")
        print(f"     ε6* range = [{vf4_eps6_star.min():.3f}, {vf4_eps6_star.max():.3f}]")
        print(f"     A[3,:] = {A[3, :]}")
        print(f"     B[3] = {B[3]}")

        return self._solve_system(A, B, "Set 2")

    def virtual_fields_set_3_debug(self):
        """Third set of virtual fields with detailed debugging"""
        X1 = self.data['X1']
        X2 = self.data['X2']
        L = self.data['L']
        w = self.data['w']
        t = self.data['t']
        F = self.data['F']
        Eps1 = self.data['Eps1']
        Eps2 = self.data['Eps2']
        Eps6 = self.data['Eps6']

        print(f"\n" + "=" * 50)
        print("VIRTUAL FIELDS SET 3")
        print("=" * 50)

        A = np.zeros((4, 4))
        B = np.zeros(4)

        # VF1: U1*=0, U2*=-x1^3 → ε1*=0, ε2*=-3x1^2, ε6*=0
        A[0, 0] = 0
        A[0, 1] = 0
        A[0, 2] = 0
        A[0, 3] = 3 * np.mean(Eps6 * X1 * X1)
        B[0] = F * L * L / (w * t)

        print(f"VF1: U1*=0, U2*=-x1³")
        print(f"     A[0,:] = {A[0, :]}")
        print(f"     B[0] = {B[0]}")

        # VF2: U1*=x1(L-x1)x2, U2*=x1^3/3-Lx1^2/2 → ε1*=(L-2x1)x2, ε2*=0, ε6*=0
        vf2_eps1_star = (L - 2 * X1) * X2
        A[1, 0] = -np.mean(Eps1 * vf2_eps1_star)
        A[1, 1] = 0
        A[1, 2] = -np.mean(Eps2 * vf2_eps1_star)
        A[1, 3] = 0
        B[1] = F * L * L / (6 * w * t)

        print(f"VF2: U1*=x1(L-x1)x2, U2*=x1³/3-Lx1²/2")
        print(f"     ε1* range = [{vf2_eps1_star.min():.3f}, {vf2_eps1_star.max():.3f}]")
        print(f"     A[1,:] = {A[1, :]}")
        print(f"     B[1] = {B[1]}")

        # VF3: U1*=0, U2*=x1(L-x1)x2 → ε1*=0, ε2*=x1(L-x1), ε6*=(L-2x1)x2
        vf3_eps2_star = X1 * (L - X1)
        vf3_eps6_star = (L - 2 * X1) * X2
        A[2, 0] = 0
        A[2, 1] = -np.mean(Eps2 * vf3_eps2_star)
        A[2, 2] = -np.mean(Eps1 * vf3_eps2_star)
        A[2, 3] = -np.mean(Eps6 * vf3_eps6_star)
        B[2] = 0

        print(f"VF3: U1*=0, U2*=x1(L-x1)x2")
        print(f"     ε2* range = [{vf3_eps2_star.min():.3f}, {vf3_eps2_star.max():.3f}]")
        print(f"     ε6* range = [{vf3_eps6_star.min():.3f}, {vf3_eps6_star.max():.3f}]")
        print(f"     A[2,:] = {A[2, :]}")
        print(f"     B[2] = {B[2]}")

        # VF4: U1*=L*sin(2πx1/L)/(2π), U2*=0 → ε1*=cos(2πx1/L), ε2*=0, ε6*=0
        vf4_eps1_star = np.cos(2 * np.pi * X1 / L)
        A[3, 0] = np.mean(Eps1 * vf4_eps1_star)
        A[3, 1] = 0
        A[3, 2] = np.mean(Eps2 * vf4_eps1_star)
        A[3, 3] = 0
        B[3] = 0

        print(f"VF4: U1*=L*sin(2πx1/L)/(2π), U2*=0")
        print(f"     ε1* range = [{vf4_eps1_star.min():.3f}, {vf4_eps1_star.max():.3f}]")
        print(f"     A[3,:] = {A[3, :]}")
        print(f"     B[3] = {B[3]}")

        return self._solve_system(A, B, "Set 3")

    def _solve_system(self, A, B, set_name):
        """Solve the linear system and display results"""
        # Check matrix condition
        cond_A = np.linalg.cond(A)
        print(f"\nMatrix condition number: {cond_A:.2e}")

        if cond_A > 1e12:
            print("WARNING: Matrix is ill-conditioned!")

        # Solve system
        try:
            Q = np.linalg.solve(A, B)

            print(f"\n{set_name} Solution:")
            print(f"Q11 = {Q[0]:.2f} MPa")
            print(f"Q22 = {Q[1]:.2f} MPa")
            print(f"Q12 = {Q[2]:.2f} MPa")
            print(f"Q66 = {Q[3]:.2f} MPa")

            # Check solution by substitution
            residual = A @ Q - B
            print(f"Residual norm: {np.linalg.norm(residual):.2e}")

            return Q, A, B

        except np.linalg.LinAlgError as e:
            print(f"Error solving system: {e}")
            return None, A, B

    def compare_with_reference(self, Q_results, Q_ref=None, material_name="reference"):
        """Compare results from all three sets with reference values"""
        if not any(Q[0] is not None for Q in Q_results):
            return

        print(f"\n" + "=" * 80)
        print("COMPARISON WITH REFERENCE VALUES")
        print("=" * 80)

        if Q_ref is None:
            print(f"No reference values provided for comparison.")
            print(f"\nComputed values (GPa):")
            print(f"{'Set':<6} {'Q11':<8} {'Q22':<8} {'Q12':<8} {'Q66':<8}")
            print("-" * 40)
            for i, (Q, _, _) in enumerate(Q_results, 1):
                if Q is not None:
                    Q_gpa = Q / 1e3
                    print(f"Set{i:<5} {Q_gpa[0]:<8.2f} {Q_gpa[1]:<8.2f} {Q_gpa[2]:<8.2f} {Q_gpa[3]:<8.2f}")
            return

        print(f"Reference material: {material_name}")
        print(f"\n{'Set':<6} {'Param':<8} {'Computed':<10} {'Reference':<10} {'Error %':<10}")
        print("-" * 50)

        param_names = ['Q11', 'Q22', 'Q12', 'Q66']

        for i, (Q, _, _) in enumerate(Q_results, 1):
            if Q is not None:
                Q_gpa = Q / 1e3
                for j, param in enumerate(param_names):
                    error_pct = abs(Q_gpa[j] - Q_ref[j]) / Q_ref[j] * 100 if Q_ref[j] != 0 else 0
                    print(f"Set{i:<5} {param:<8} {Q_gpa[j]:<10.2f} {Q_ref[j]:<10.2f} {error_pct:<10.3f}")
                print("-" * 50)

        # Summary comparison table
        print(f"\nSUMMARY COMPARISON TABLE")
        print(f"{'Parameter':<8} {'Set 1':<10} {'Set 2':<10} {'Set 3':<10} {'Reference':<10}")
        print("-" * 60)

        for j, param in enumerate(param_names):
            print(f"{param:<8}", end="")
            for i, (Q, _, _) in enumerate(Q_results):
                if Q is not None:
                    Q_gpa = Q / 1e3
                    print(f" {Q_gpa[j]:<10.2f}", end="")
                else:
                    print(f" {'N/A':<10}", end="")
            if Q_ref is not None:
                print(f" {Q_ref[j]:<10.2f}")
            else:
                print()

    def run_analysis(self, Q_ref=None, material_name="reference"):
        """Run debug analysis with all three virtual field sets"""
        print("=" * 80)
        print("VFM DEBUG ANALYSIS - THREE VIRTUAL FIELD SETS")
        print("=" * 80)

        # Check virtual field consistency
        self.check_virtual_field_consistency()

        # Run all three sets
        Q1, A1, B1 = self.virtual_fields_set_1_debug()
        Q2, A2, B2 = self.virtual_fields_set_2_debug()
        Q3, A3, B3 = self.virtual_fields_set_3_debug()

        # Store results
        Q_results = [(Q1, A1, B1), (Q2, A2, B2), (Q3, A3, B3)]

        # Compare with reference
        self.compare_with_reference(Q_results, Q_ref, material_name)

        return Q_results


    def plot_strain_fields(self, save_path='strain_fields.png'):
        """Plot the strain field distributions with enhanced formatting"""

        FS = 22
        cor = 'BrBG'
        # Enable LaTeX rendering
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.size'] = FS  # Base font size

        # Get data
        X1 = self.data['X1']
        X2 = self.data['X2']
        strain_data = [self.data['Eps1'], self.data['Eps2'], self.data['Eps6']]
        strain_labels = [r'$\varepsilon_1$', r'$\varepsilon_2$', r'$\varepsilon_6$']

        # Create figure
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))

        # Plot each strain component
        for i, (data, label, ax) in enumerate(zip(strain_data, strain_labels, axes)):
            # Use actual data range for colorbar
            vmin, vmax = np.min(data), np.max(data)

            # Create contour plot
            im = ax.contourf(X1, X2, data, levels=20, cmap=cor,
                             vmin=vmin, vmax=vmax)
            # correct aspect ratio
            ax.set_aspect('equal')
            # Formatting
            ax.set_title(label, fontsize=FS + 20)
            ax.set_xlabel(r'$x_1$ (mm)', fontsize=FS)
            ax.set_ylabel(r'$x_2$ (mm)', fontsize=FS)
            ax.tick_params(labelsize=FS - 4)

            # Rectangular colorbar
            cbar = plt.colorbar(im, ax=ax, shrink=0.6, aspect=30, pad=0.05, fraction=0.04)
            cbar.ax.tick_params(labelsize=FS - 4)
            cbar.formatter.set_powerlimits((0, 0))
            cbar.update_ticks()


        plt.tight_layout()

        # Save the figure with high DPI
        plt.savefig(save_path, dpi=150, bbox_inches='tight',
                    facecolor='white', edgecolor='none')

        plt.show()

def main():
    """Main function - Edit reference values directly in this section"""

    # ========== EDIT REFERENCE VALUES HERE ==========
    # Change these values for your specific material:
    Q_ref = np.array([15.536, 1.965, 0.926, 1.109])  # Q11, Q22, Q12, Q66 in GPa (from your figure)
    material_name = "Wood"

    # Set to None if you don't want comparison:
    # Q_ref = None
    # material_name = "No reference"
    # ===============================================

    vfm = IosipescuVFM(data_source='csv')
    vfm.load_data('.')  # Load from current directory

    Q_results = vfm.run_analysis(Q_ref, material_name)

    # Plot strain fields
    try:
        vfm.plot_strain_fields()
    except Exception as e:
        print(f"Could not plot strain fields: {e}")

if __name__ == "__main__":
    main()