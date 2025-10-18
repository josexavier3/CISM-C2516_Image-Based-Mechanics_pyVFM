import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


class VirtualFieldsMethod:
    """
    Virtual Fields Method implementation for disc under compression
    Includes Set 1, Set 2, Set 3 virtual fields and combined least-squares solution
    WITH NOISE SENSITIVITY ANALYSIS
    """

    def __init__(self, data_file, force_magnitude, disc_diameter, thickness):
        """
        Initialize VFM analysis

        Parameters:
        - data_file: path to ANSYS results file
        - force_magnitude: applied force F (N)
        - disc_diameter: disc diameter d (mm)
        - thickness: disc thickness t (mm)
        """
        self.F = force_magnitude
        self.d = disc_diameter
        self.t = thickness
        self.data = self.load_ansys_data(data_file)
        self.data_clean = self.data.copy() if self.data is not None else None  # Store clean data

    def load_ansys_data(self, filename):
        """
        Load ANSYS strain data with area information (CSV format, semicolon-separated)
        Parse CSV format: Area_mm2; X_Coord; Y_Coord; U_X; U_Y; Strain_X; Strain_Y; Strain_XY
        """
        try:
            # Read the raw file to understand the format
            with open(filename, 'r') as f:
                lines = f.readlines()

            # Find the data start (skip header lines)
            data_start = 0
            for i, line in enumerate(lines):
                # Look for Area_mm2 in header (semicolon-separated)
                if 'Area_mm2' in line and 'X_Coord' in line:
                    data_start = i + 1  # Skip header line
                    break

            print(f"Found data starting at line {data_start}")
            if data_start > 0:
                print(f"Header line: {lines[data_start - 1].strip()}")

            # Parse data lines
            parsed_data = []
            for line_num, line in enumerate(lines[data_start:], start=data_start):
                line = line.strip()
                if not line:
                    continue

                # Split by semicolon and clean up whitespace
                parts = line.split(';')
                # Strip whitespace from each part
                parts = [p.strip() for p in parts]

                if len(parts) < 8:  # Need at least 8 columns: Area, X, Y, UX, UY, StrainX, StrainY, StrainXY
                    continue

                try:
                    # Format: Area_mm2, X_Coord, Y_Coord, U_X, U_Y, Strain_X, Strain_Y, Strain_XY
                    area = self.parse_ansys_number(parts[0])
                    x = float(parts[1])
                    y = float(parts[2])

                    # Strains (columns 5, 6, 7)
                    exx = self.parse_ansys_number(parts[5])  # Strain_X
                    eyy = self.parse_ansys_number(parts[6])  # Strain_Y
                    exy = self.parse_ansys_number(parts[7])  # Strain_XY (engineering shear strain)

                    parsed_data.append([x, y, exx, eyy, exy, area])

                except (ValueError, IndexError) as e:
                    if line_num < data_start + 5:  # Show first few parsing errors
                        print(f"Warning: Could not parse line {line_num}: {line[:50]}...")
                    continue

            if not parsed_data:
                print("No valid data found. Please check the file format.")
                print("First few lines of file:")
                for i, line in enumerate(lines[:10]):
                    print(f"{i}: {line.strip()}")
                return None

            # Create DataFrame
            data = pd.DataFrame(parsed_data, columns=['x', 'y', 'exx', 'eyy', 'exy', 'area'])

            print(f"Successfully loaded {len(data)} data points")
            print(f"Data columns: {data.columns.tolist()}")
            print(f"Data preview:\n{data.head()}")
            print(f"Data stats:\n{data.describe()}")
            print(f"Total area from elements: {data['area'].sum():.2f} mm²")

            return data

        except Exception as e:
            print(f"Error loading data: {e}")
            print("Please check your file format.")
            return None

    def parse_ansys_number(self, num_str):
        """
        Parse ANSYS scientific notation (e.g., '0.3918E-04' -> 3.918e-5)
        """
        try:
            # Handle ANSYS E notation by replacing with standard e notation
            num_str = num_str.replace('E', 'e').replace('D', 'e')
            return float(num_str)
        except:
            return 0.0

    def add_gaussian_noise(self, noise_level):
        """
        Add Gaussian white noise to strain measurements

        Parameters:
        - noise_level: amplitude of noise (e.g., 1e-6, 1e-5, 1e-4)
        """
        if self.data_clean is None:
            print("No clean data available")
            return

        # Reset to clean data
        self.data = self.data_clean.copy()

        # Add Gaussian noise to strain components
        np.random.seed(42)  # For reproducibility
        n_points = len(self.data)

        noise_exx = np.random.normal(0, noise_level, n_points)
        noise_eyy = np.random.normal(0, noise_level, n_points)
        noise_exy = np.random.normal(0, noise_level, n_points)

        self.data['exx'] = self.data_clean['exx'] + noise_exx
        self.data['eyy'] = self.data_clean['eyy'] + noise_eyy
        self.data['exy'] = self.data_clean['exy'] + noise_exy

        print(f"Added Gaussian noise with amplitude {noise_level:.1e}")

    def compute_virtual_fields_set1(self):
        """
        SET 1: Simple virtual fields
        VF1: Pure compression (u1*=0, u2*=-k1*x2)
        VF2: Pure horizontal extension (u1*=k2*x1, u2*=0)
        """
        x = self.data['x'].values
        y = self.data['y'].values

        # Virtual Field 1: Pure compression
        k1 = 1.0
        vf1_strains = {
            'e1_star': np.zeros(len(x)),
            'e2_star': -k1 * np.ones(len(x)),
            'e6_star': np.zeros(len(x))
        }

        # Virtual Field 2: Pure horizontal extension
        k2 = 1.0
        vf2_strains = {
            'e1_star': k2 * np.ones(len(x)),
            'e2_star': np.zeros(len(x)),
            'e6_star': np.zeros(len(x))
        }

        return [vf1_strains, vf2_strains], [k1, k2]

    def compute_virtual_fields_set2(self):
        """
        SET 2: Trigonometric and Exponential virtual fields
        VF3: Trigonometric (u1*=0, u2*=-sin(πx2/2d))
        VF4: Exponential (u1*=exp(x1/d), u2*=0)
        """
        x = self.data['x'].values
        y = self.data['y'].values

        # Virtual Field 3: Trigonometric
        vf3_strains = {
            'e1_star': np.zeros(len(x)),
            'e2_star': -(np.pi / (2 * self.d)) * np.cos(np.pi * y / (2 * self.d)),
            'e6_star': np.zeros(len(x))
        }

        # Virtual Field 4: Exponential
        vf4_strains = {
            'e1_star': (1 / self.d) * np.exp(x / self.d),
            'e2_star': np.zeros(len(x)),
            'e6_star': np.zeros(len(x))
        }

        return [vf3_strains, vf4_strains]

    def compute_virtual_fields_set3(self):
        """
        SET 3: Higher-order polynomial and exponential
        VF5: Cubic polynomial (u1*=0, u2*=-x2^3)
        VF6: Quadratic exponential (u1*=exp(x1^2/d), u2*=0)
        """
        x = self.data['x'].values
        y = self.data['y'].values

        # Virtual Field 5: Higher-order polynomial
        # e1* = 0, e2* = -3x2^2, e6* = 0
        vf5_strains = {
            'e1_star': np.zeros(len(x)),
            'e2_star': -3 * y ** 2,
            'e6_star': np.zeros(len(x))
        }

        # Virtual Field 6: Higher-order exponential
        # e1* = (2x1/d)exp(x1^2/d), e2* = 0, e6* = 0
        vf6_strains = {
            'e1_star': (2 * x / self.d) * np.exp(x ** 2 / self.d),
            'e2_star': np.zeros(len(x)),
            'e6_star': np.zeros(len(x))
        }

        return [vf5_strains, vf6_strains]

    def compute_integrals_general(self, vf_strains_list):
        """
        Compute the A matrix for any set of virtual fields

        Parameters:
        - vf_strains_list: list of virtual field strain dictionaries

        Returns:
        - A: coefficient matrix (n_vf x 2) where n_vf is number of virtual fields
        """
        # Extract measured strains
        e1 = self.data['exx'].values  # strain in x-direction
        e2 = self.data['eyy'].values  # strain in y-direction

        # Get actual element areas from the data
        Sa = self.data['area'].values  # Area for each element in mm²

        n_vf = len(vf_strains_list)
        A = np.zeros((n_vf, 2))

        for i, vf_strains in enumerate(vf_strains_list):
            A[i, 0] = np.sum(Sa * (e2 * vf_strains['e2_star'] + e1 * vf_strains['e1_star']))
            A[i, 1] = np.sum(Sa * (e1 * vf_strains['e2_star'] + e2 * vf_strains['e1_star']))

        return A

    def compute_external_work_set1(self, k1, k2):
        """Compute external virtual work for Set 1 virtual fields"""
        B1 = self.F * k1 * self.d / self.t
        B2 = 0.0
        return np.array([B1, B2])

    def compute_external_work_set2(self):
        """Compute external virtual work for Set 2 virtual fields"""
        B3 = -self.F * (-1) / self.t
        B4 = 0.0
        return np.array([B3, B4])

    def compute_external_work_set3(self):
        """Compute external virtual work for Set 3 virtual fields"""
        B5 = -self.F * (-self.d ** 3) / self.t
        B6 = 0.0
        return np.array([B5, B6])

    def solve_system(self, A, B, method='direct'):
        """
        Solve the linear system AQ = B for material parameters Q11, Q12

        Parameters:
        - A: coefficient matrix
        - B: external work vector
        - method: 'direct' for square system, 'least_squares' for overdetermined

        Returns:
        - Q11, Q12: material parameters
        - cond: condition number of matrix A
        """
        if method == 'direct':
            try:
                # Compute condition number
                cond = np.linalg.cond(A)

                Q = np.linalg.solve(A, B)
                Q11, Q12 = Q[0], Q[1]
                return Q11, Q12, cond
            except np.linalg.LinAlgError:
                print("Error: Singular matrix - cannot solve system")
                return None, None, np.inf

        elif method == 'least_squares':
            # Compute condition number
            cond = np.linalg.cond(A)

            Q = np.linalg.lstsq(A, B, rcond=None)[0]
            Q11, Q12 = Q[0], Q[1]
            return Q11, Q12, cond

    def convert_to_engineering_constants(self, Q11, Q12, verbose=False):
        """Convert stiffness matrix components to engineering constants"""
        if Q11 is None or Q12 is None:
            return None, None

        nu = Q12 / Q11
        E = Q11 * (1 - nu ** 2)

        if verbose:
            print(f"\nEngineering constants:")
            print(f"  Young's modulus E = {E:.2f} MPa")
            print(f"  Poisson's ratio ν = {nu:.4f}")

        return E, nu

    def run_analysis(self, E_ref=None, nu_ref=None):
        """Run complete VFM analysis with optional reference values"""
        if self.data is None:
            return

        print("\n" + "=" * 80)
        print("VIRTUAL FIELDS METHOD ANALYSIS")
        print("=" * 80)
        results = {}

        # SET 1
        print("\nSET 1: Simple Virtual Fields")
        vf1, k1 = self.compute_virtual_fields_set1()
        A1 = self.compute_integrals_general(vf1)
        B1 = self.compute_external_work_set1(k1[0], k1[1])
        Q11_s1, Q12_s1, cond1 = self.solve_system(A1, B1)
        if Q11_s1:
            E_s1, nu_s1 = self.convert_to_engineering_constants(Q11_s1, Q12_s1, verbose=True)
            results['Set1'] = {'E': E_s1, 'nu': nu_s1, 'Q11': Q11_s1, 'Q12': Q12_s1, 'cond': cond1}

        # SET 2
        print("\nSET 2: Trigonometric & Exponential")
        vf2 = self.compute_virtual_fields_set2()
        A2 = self.compute_integrals_general(vf2)
        B2 = self.compute_external_work_set2()
        Q11_s2, Q12_s2, cond2 = self.solve_system(A2, B2)
        if Q11_s2:
            E_s2, nu_s2 = self.convert_to_engineering_constants(Q11_s2, Q12_s2, verbose=True)
            results['Set2'] = {'E': E_s2, 'nu': nu_s2, 'Q11': Q11_s2, 'Q12': Q12_s2, 'cond': cond2}

        # SET 3
        print("\nSET 3: Higher-order Polynomial & Exponential")
        vf3 = self.compute_virtual_fields_set3()
        A3 = self.compute_integrals_general(vf3)
        B3 = self.compute_external_work_set3()
        Q11_s3, Q12_s3, cond3 = self.solve_system(A3, B3)
        if Q11_s3:
            E_s3, nu_s3 = self.convert_to_engineering_constants(Q11_s3, Q12_s3, verbose=True)
            results['Set3'] = {'E': E_s3, 'nu': nu_s3, 'Q11': Q11_s3, 'Q12': Q12_s3, 'cond': cond3}

        # COMBINED
        print("\nCOMBINED: Least-Squares")
        vf_all = vf1 + vf2 + vf3
        A_all = self.compute_integrals_general(vf_all)
        B_all = np.concatenate([B1, B2, B3])
        Q11_c, Q12_c, cond_c = self.solve_system(A_all, B_all, method='least_squares')
        if Q11_c:
            E_c, nu_c = self.convert_to_engineering_constants(Q11_c, Q12_c, verbose=True)
            results['Combined'] = {'E': E_c, 'nu': nu_c, 'Q11': Q11_c, 'Q12': Q12_c, 'cond': cond_c}

        # SUMMARY
        print("\n" + "=" * 80)
        print("SUMMARY")
        print("=" * 80)

        if E_ref is not None and nu_ref is not None:
            print(f"Reference values: E = {E_ref:.2f} GPa, ν = {nu_ref:.4f}\n")
            print(f"{'Method':<20} {'E (GPa)':<25} {'ν':<25} {'Condition #':<15}")
            print("-" * 85)
            for key, v in results.items():
                E_GPa = v['E'] / 1000.0
                E_error = ((E_GPa - E_ref) / E_ref) * 100
                nu_error = ((v['nu'] - nu_ref) / nu_ref) * 100
                E_str = f"{E_GPa:.4f} ({E_error:+.2f}%)"
                nu_str = f"{v['nu']:.6f} ({nu_error:+.2f}%)"
                cond_str = f"{v['cond']:.2e}"
                print(f"{key:<20} {E_str:<25} {nu_str:<25} {cond_str:<15}")
        else:
            print(f"{'Method':<20} {'E (GPa)':<12} {'ν':<12} {'Condition #':<15}")
            print("-" * 60)
            for key, v in results.items():
                E_GPa = v['E'] / 1000.0
                cond_str = f"{v['cond']:.2e}"
                print(f"{key:<20} {E_GPa:<12.4f} {v['nu']:<12.6f} {cond_str:<15}")

        return results

    def noise_sensitivity_analysis(self, noise_levels, E_ref=210.0, nu_ref=0.3):
        """
        Perform noise sensitivity analysis across multiple noise levels

        Parameters:
        - noise_levels: list of noise amplitudes (e.g., [1e-6, 1e-5, 1e-4])
        - E_ref: reference Young's modulus in GPa
        - nu_ref: reference Poisson's ratio
        """
        print("\n" + "=" * 100)
        print("NOISE SENSITIVITY & CONDITIONING ANALYSIS")
        print("=" * 100)

        results_table = []

        for noise_level in noise_levels:
            print(f"\n{'=' * 100}")
            print(f"NOISE LEVEL: {noise_level:.0e}")
            print(f"{'=' * 100}")

            # Add noise to data
            self.add_gaussian_noise(noise_level)

            # Run analysis
            results = {}

            # SET 1
            vf1, k1 = self.compute_virtual_fields_set1()
            A1 = self.compute_integrals_general(vf1)
            B1 = self.compute_external_work_set1(k1[0], k1[1])
            Q11_s1, Q12_s1, cond1 = self.solve_system(A1, B1)
            if Q11_s1:
                E_s1, nu_s1 = self.convert_to_engineering_constants(Q11_s1, Q12_s1)
                results['Set1'] = {'E': E_s1 / 1000, 'nu': nu_s1, 'cond': cond1}

            # SET 2
            vf2 = self.compute_virtual_fields_set2()
            A2 = self.compute_integrals_general(vf2)
            B2 = self.compute_external_work_set2()
            Q11_s2, Q12_s2, cond2 = self.solve_system(A2, B2)
            if Q11_s2:
                E_s2, nu_s2 = self.convert_to_engineering_constants(Q11_s2, Q12_s2)
                results['Set2'] = {'E': E_s2 / 1000, 'nu': nu_s2, 'cond': cond2}

            # SET 3
            vf3 = self.compute_virtual_fields_set3()
            A3 = self.compute_integrals_general(vf3)
            B3 = self.compute_external_work_set3()
            Q11_s3, Q12_s3, cond3 = self.solve_system(A3, B3)
            if Q11_s3:
                E_s3, nu_s3 = self.convert_to_engineering_constants(Q11_s3, Q12_s3)
                results['Set3'] = {'E': E_s3 / 1000, 'nu': nu_s3, 'cond': cond3}

            # Store results
            results_table.append({
                'noise_level': noise_level,
                'Set1_E': results.get('Set1', {}).get('E', np.nan),
                'Set1_nu': results.get('Set1', {}).get('nu', np.nan),
                'Set1_cond': results.get('Set1', {}).get('cond', np.nan),
                'Set2_E': results.get('Set2', {}).get('E', np.nan),
                'Set2_nu': results.get('Set2', {}).get('nu', np.nan),
                'Set2_cond': results.get('Set2', {}).get('cond', np.nan),
                'Set3_E': results.get('Set3', {}).get('E', np.nan),
                'Set3_nu': results.get('Set3', {}).get('nu', np.nan),
                'Set3_cond': results.get('Set3', {}).get('cond', np.nan),
            })

        # Display summary table
        print("\n" + "=" * 100)
        print("NOISE SENSITIVITY SUMMARY")
        print("=" * 100)
        print(f"\n{'Metric':<15} {'Level':<12} {'Set 1':<20} {'Set 2':<20} {'Set 3':<20}")
        print("-" * 100)

        for row in results_table:
            nl = row['noise_level']
            print(f"{'E (GPa)':<15} {nl:<12.0e} {row['Set1_E']:>20.3f} {row['Set2_E']:>20.3f} {row['Set3_E']:>20.3f}")

        print()
        for row in results_table:
            nl = row['noise_level']
            print(f"{'ν':<15} {nl:<12.0e} {row['Set1_nu']:>20.4f} {row['Set2_nu']:>20.4f} {row['Set3_nu']:>20.4f}")

        print()
        for row in results_table:
            nl = row['noise_level']
            print(
                f"{'Condition #':<15} {nl:<12.0e} {row['Set1_cond']:>20.2e} {row['Set2_cond']:>20.2e} {row['Set3_cond']:>20.2e}")

        # Calculate errors
        print("\n" + "=" * 100)
        print("ERRORS RELATIVE TO REFERENCE VALUES")
        print("=" * 100)
        print(f"Reference: E = {E_ref} GPa, ν = {nu_ref}\n")

        print(f"{'Metric':<15} {'Level':<12} {'Set 1 Error (%)':<20} {'Set 2 Error (%)':<20} {'Set 3 Error (%)':<20}")
        print("-" * 100)

        for row in results_table:
            nl = row['noise_level']
            E1_err = ((row['Set1_E'] - E_ref) / E_ref) * 100
            E2_err = ((row['Set2_E'] - E_ref) / E_ref) * 100
            E3_err = ((row['Set3_E'] - E_ref) / E_ref) * 100
            print(f"{'E error':<15} {nl:<12.0e} {E1_err:>+20.2f} {E2_err:>+20.2f} {E3_err:>+20.2f}")

        print()
        for row in results_table:
            nl = row['noise_level']
            nu1_err = ((row['Set1_nu'] - nu_ref) / nu_ref) * 100
            nu2_err = ((row['Set2_nu'] - nu_ref) / nu_ref) * 100
            nu3_err = ((row['Set3_nu'] - nu_ref) / nu_ref) * 100
            print(f"{'ν error':<15} {nl:<12.0e} {nu1_err:>+20.2f} {nu2_err:>+20.2f} {nu3_err:>+20.2f}")

        return results_table

    def plot_strain_fields(self, save_figure=True, filename='strain_fields.png', dpi=300, noise_level=None):
        """
        Plot strain fields with optional noise

        Parameters:
        - save_figure: bool, whether to save the figure
        - filename: str, output filename
        - dpi: int, resolution
        - noise_level: float or None, if provided, adds Gaussian noise before plotting
        """
        if self.data is None:
            return

        # Add noise if requested
        if noise_level is not None:
            self.add_gaussian_noise(noise_level)
            print(f"Plotting strain fields with noise level: {noise_level:.0e}")
        else:
            # Reset to clean data if no noise requested
            if self.data_clean is not None:
                self.data = self.data_clean.copy()
            print("Plotting clean strain fields (no noise)")

        # Enable LaTeX
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        FS = 18
        cor = 'BrBG'

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), facecolor='#fdf6e3')

        x = self.data['x'].values
        y = self.data['y'].values
        exx = self.data['exx'].values
        eyy = self.data['eyy'].values

        # Grid
        xi = np.linspace(x.min(), x.max(), 200)
        yi = np.linspace(y.min(), y.max(), 200)
        Xi, Yi = np.meshgrid(xi, yi)

        # Interpolate
        exx_i = griddata((x, y), exx, (Xi, Yi), method='cubic')
        eyy_i = griddata((x, y), eyy, (Xi, Yi), method='cubic')

        # Mask
        R = self.d / 2
        mask = Xi ** 2 + (Yi - self.d / 2) ** 2 > R ** 2
        exx_i[mask] = np.nan
        eyy_i[mask] = np.nan

        # Add noise level to title if applicable
        title_suffix = f' (noise: {noise_level:.0e})' if noise_level is not None else ''

        # Plot exx
        c1 = ax1.contourf(Xi, Yi, exx_i, levels=15, cmap=cor)
        ax1.contour(Xi, Yi, exx_i, levels=15, colors='k', linewidths=0.3, alpha=0.3)
        ax1.set_title(r'Strain $\varepsilon_{xx}$' + title_suffix, fontsize=FS + 8)
        ax1.set_xlabel(r'$x$ (mm)', fontsize=FS)
        ax1.set_ylabel(r'$y$ (mm)', fontsize=FS)
        ax1.tick_params(labelsize=FS)
        ax1.set_aspect('equal')
        cbar1 = plt.colorbar(c1, ax=ax1, fraction=0.046, pad=0.04)
        cbar1.ax.tick_params(labelsize=FS + 2)

        # Plot eyy
        c2 = ax2.contourf(Xi, Yi, eyy_i, levels=15, cmap=cor)
        ax2.contour(Xi, Yi, eyy_i, levels=15, colors='k', linewidths=0.3, alpha=0.3)
        ax2.set_title(r'Strain $\varepsilon_{yy}$' + title_suffix, fontsize=FS + 8)
        ax2.set_xlabel(r'$x$ (mm)', fontsize=FS)
        ax2.set_ylabel(r'$y$ (mm)', fontsize=FS)
        ax2.tick_params(labelsize=FS)
        ax2.set_aspect('equal')
        cbar2 = plt.colorbar(c2, ax=ax2, fraction=0.046, pad=0.04)
        cbar2.ax.tick_params(labelsize=FS + 2)

        plt.tight_layout()

        if save_figure:
            fig.savefig(filename, dpi=dpi, bbox_inches='tight')
            print(f"Saved: {filename}")

        plt.rcParams['text.usetex'] = False
        plt.show()
        return fig

# Example usage
if __name__ == "__main__":
    # Define test parameters
    force = 47981.1  # Reaction force in N
    diameter = 200.0  # Disc diameter in mm
    thickness = 3.0  # Disc thickness in mm
    E_ref = 210.0  # Reference Young's modulus in GPa
    nu_ref = 0.3  # Reference Poisson's ratio

    # Initialize VFM analysis
    vfm = VirtualFieldsMethod(
        data_file='FEM2VFM.csv',
        force_magnitude=force,
        disc_diameter=diameter,
        thickness=thickness
    )

    try:
        # Plot strain fields
        vfm.plot_strain_fields(save_figure=True, filename='disc_strain_fields_noise.png', dpi=150)

        vfm.plot_strain_fields(
            save_figure=True,
            filename='disc_strain_fields_noise_1e-6.png',
            dpi=150,
            noise_level=1e-6
        )

        vfm.plot_strain_fields(
            save_figure=True,
            filename='disc_strain_fields_noise_1e-5.png',
            dpi=150,
            noise_level=1e-5
        )

        # Plot strain fields with noise level 1e-4
        vfm.plot_strain_fields(
            save_figure=True,
            filename='disc_strain_fields_noise_1e-4.png',
            dpi=150,
            noise_level=1e-4
        )

        # Run analysis without noise (clean data)
        print("\n" + "=" * 100)
        print("ANALYSIS WITH CLEAN DATA (NO NOISE)")
        print("=" * 100)
        results_clean = vfm.run_analysis(E_ref=E_ref, nu_ref=nu_ref)

        # Noise sensitivity analysis
        noise_levels = [1e-6, 1e-5, 1e-4]
        results_noise = vfm.noise_sensitivity_analysis(noise_levels, E_ref=E_ref, nu_ref=nu_ref)

    except Exception as e:
        print(f"Analysis failed: {e}")
        import traceback

        traceback.print_exc()