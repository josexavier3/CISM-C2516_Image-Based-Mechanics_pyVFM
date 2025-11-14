"""
Course reference:
CISM Advanced School: "Image-Based Mechanics: An Overview of Experimental
and Numerical Approaches"
Udine, Italy, October 6-10, 2025
Coordinated by: Julien Réthoré and José Xavier

Lecture topic:
"2D & Stereo Digital Image Correlation: Guidance and Practical Concepts"

Developed by:
José Xavier
Universidade NOVA de Lisboa, NOVA FCT, UNIDEMI
https://userweb.fct.unl.pt/~jmc.xavier/index.html

For more information: https://cism.it/en/activities/courses/C2516/

=============================================================================

DIC Speckle Pattern Calculator (Pitch & Resolution Driven)
----------------------------------------------------------
User provides ONE of the optical parameters (f, s, or FOV_width).
FOV_height is AUTOMATICALLY calculated from sensor aspect ratio.

Key symbols (units in brackets):
  - res_w_px, res_h_px  : image resolution along width/height [px]
  - pixel_pitch_um      : pixel pitch on sensor [µm]
  - sensor_w_mm, sensor_h_mm : sensor physical size [mm]
  - f                   : lens focal length [mm]
  - s                   : object (working) distance: object plane → lens principal plane [mm]
  - m                   : magnification (image size / object size) [-]
  - FOV_w, FOV_h        : field-of-view along width/height on object [mm]
  - p_obj               : pixel size on the object [mm]
  - scale               : image scale [px/mm]

Core equations:
  (1) sensor_w_mm  = res_w_px * pixel_pitch_um / 1000
  (2) m = f / (s - f)
  (3) FOV_w = sensor_w_mm / m; FOV_h = sensor_h_mm / m
  (4) m = sensor_w_mm / FOV_w
  (5) s = f * (1 + 1/m)
  (6) f = (m * s) / (1 + m)
  (7) p_obj = pixel_pitch_um / (m * 1000)
  (8) scale = m * 1000 / pixel_pitch_um
  (9) speckle_size_mm = N / scale

Working Distance Calculator - Multiple Optical Models
======================================================

Computes working distance (object distance) using various optical models
when only focal length is provided. Useful for comparing against commercial
DIC software that may use non-standard definitions.

Since lens specifications are rarely available, this module provides
multiple models with typical correction factors based on common lens
designs and mounting standards.

Reference Distance Definitions:
  s_principal  : Distance from object to lens PRINCIPAL PLANE (thin-lens model)
  s_front      : Distance from object to FRONT LENS SURFACE
  s_BFD        : Distance from object to sensor, accounting for BACK FOCAL DISTANCE
  s_flange     : Distance from object to CAMERA MOUNT REFERENCE PLANE (flange)
"""

"""
Course reference:
CISM Advanced School: "Image-Based Mechanics: An Overview of Experimental
and Numerical Approaches"
Udine, Italy, October 6-10, 2025
Coordinated by: Julien Réthoré and José Xavier

Lecture topic:
"2D & Stereo Digital Image Correlation: Guidance and Practical Concepts"

Developed by:
José Xavier
Universidade NOVA de Lisboa, NOVA FCT, UNIDEMI
https://userweb.fct.unl.pt/~jmc.xavier/index.html

For more information: https://cism.it/en/activities/courses/C2516/

=============================================================================

DIC Speckle Pattern Calculator (Pitch & Resolution Driven)
----------------------------------------------------------
User provides ONE of the optical parameters (f, s, or FOV_width).
FOV_height is AUTOMATICALLY calculated from sensor aspect ratio.

Key symbols (units in brackets):
  - res_w_px, res_h_px  : image resolution along width/height [px]
  - pixel_pitch_um      : pixel pitch on sensor [µm]
  - sensor_w_mm, sensor_h_mm : sensor physical size [mm]
  - f                   : lens focal length [mm]
  - s                   : object (working) distance: object plane → lens principal plane [mm]
  - m                   : magnification (image size / object size) [-]
  - FOV_w, FOV_h        : field-of-view along width/height on object [mm]
  - p_obj               : pixel size on the object [mm]
  - scale               : image scale [px/mm]

Core equations:
  (1) sensor_w_mm  = res_w_px * pixel_pitch_um / 1000
  (2) m = f / (s - f)
  (3) FOV_w = sensor_w_mm / m; FOV_h = sensor_h_mm / m
  (4) m = sensor_w_mm / FOV_w
  (5) s = f * (1 + 1/m)
  (6) f = (m * s) / (1 + m)
  (7) p_obj = pixel_pitch_um / (m * 1000)
  (8) scale = m * 1000 / pixel_pitch_um
  (9) speckle_size_mm = N / scale
"""

import numpy as np

# ==============================================================================
# PAINTING TECHNIQUES (speckle size windows in mm)
# ==============================================================================
PAINTING_TECHNIQUES = {
    'Fine airbrush (0.3mm nozzle)': (0.05, 0.30),
    'Standard airbrush (0.5mm)': (0.20, 0.80),
    'Spray paint (fine)': (0.40, 1.50),
    'Spray paint (coarse)': (0.80, 3.00),
    'Speckle stamp': (0.30, 1.50),
    'Toner powder': (0.01, 0.10),
    'Lithography/etching': (0.001, 0.050),
}


# ==============================================================================
# BASIC OPTICS HELPERS
# ==============================================================================

def calculate_magnification(working_distance_mm, focal_length_mm):
    """Eq. (2) m = f / (s - f)"""
    return focal_length_mm / (working_distance_mm - focal_length_mm)


def calculate_field_of_view(sensor_dim_mm, magnification):
    """Eq. (3) FOV_dim = sensor_dim / m"""
    return sensor_dim_mm / magnification


def calculate_pixel_size_on_object(pixel_pitch_um, magnification):
    """Eq. (7) p_obj = pixel_pitch_um / (m * 1000)"""
    return pixel_pitch_um / (magnification * 1000.0)


def calculate_image_scale(magnification, pixel_pitch_um):
    """Eq. (8) scale = m * 1000 / pixel_pitch_um"""
    return magnification * 1000.0 / pixel_pitch_um


def calculate_speckle_size(pixels_per_speckle, image_scale_px_per_mm):
    """Eq. (9) speckle_size_mm = N / scale"""
    return pixels_per_speckle / image_scale_px_per_mm


def recommend_technique(speckle_size_mm):
    rec = []
    for technique, (mn, mx) in PAINTING_TECHNIQUES.items():
        if mn <= speckle_size_mm <= mx:
            rec.append(technique)
    return rec if rec else ["Custom technique needed"]


# ==============================================================================
# FLEXIBLE OPTICS SOLVER
# ==============================================================================

def solve_optics_with_auto_fov(sensor_w_mm_native, sensor_h_mm_native,
                               focal_length_mm=None,
                               working_distance_mm=None,
                               fov_width_mm=None):
    """
    Solve for optical quantities. User provides EXACTLY TWO of:
      - focal_length_mm (f)
      - working_distance_mm (s)
      - fov_width_mm (FOV along longest sensor side)

    FOV_height is AUTOMATICALLY calculated from sensor aspect ratio:
      fov_height_mm = fov_width_mm * (sensor_h_short / sensor_w_long)

    Returns:
      f, s, FOV_w (long side), FOV_h (short side, auto-calculated),
      and the ordered sensor sizes (sensor_w_mm_long, sensor_h_mm_short)
    """

    # Determine which sensor side is longer (define "width" for FOV convention)
    if sensor_w_mm_native >= sensor_h_mm_native:
        sensor_w_mm_long = sensor_w_mm_native
        sensor_h_mm_short = sensor_h_mm_native
    else:
        sensor_w_mm_long = sensor_h_mm_native
        sensor_h_mm_short = sensor_w_mm_native

    provided = {
        'f': focal_length_mm is not None,
        's': working_distance_mm is not None,
        'FOV_w': fov_width_mm is not None,
    }

    if sum(provided.values()) != 2:
        raise ValueError("Provide exactly TWO of: focal_length_mm, working_distance_mm, fov_width_mm")

    # ----------------------------
    # Case 1: f & s -> FOV via m
    # ----------------------------
    if (focal_length_mm is not None) and (working_distance_mm is not None):
        if working_distance_mm <= focal_length_mm:
            raise ValueError("Working distance must be greater than focal length.")
        m = calculate_magnification(working_distance_mm, focal_length_mm)
        fov_width_mm = calculate_field_of_view(sensor_w_mm_long, m)
        # Auto-calculate FOV_h from aspect ratio
        fov_height_mm = fov_width_mm * (sensor_h_mm_short / sensor_w_mm_long)

    # ----------------------------
    # Case 2: f & FOV_w -> s
    # ----------------------------
    elif (focal_length_mm is not None) and (fov_width_mm is not None):
        m = sensor_w_mm_long / fov_width_mm
        working_distance_mm = focal_length_mm * (1.0 + 1.0 / m)
        # Auto-calculate FOV_h from aspect ratio
        fov_height_mm = fov_width_mm * (sensor_h_mm_short / sensor_w_mm_long)

    # ----------------------------
    # Case 3: s & FOV_w -> f
    # ----------------------------
    elif (working_distance_mm is not None) and (fov_width_mm is not None):
        m = sensor_w_mm_long / fov_width_mm
        focal_length_mm = (m * working_distance_mm) / (1.0 + m)
        # Auto-calculate FOV_h from aspect ratio
        fov_height_mm = fov_width_mm * (sensor_h_mm_short / sensor_w_mm_long)

    else:
        raise ValueError("Unsupported combination. Provide exactly two valid quantities.")

    # Physical feasibility check
    if working_distance_mm <= focal_length_mm:
        raise ValueError("Computed working distance ≤ focal length (non-physical). Check inputs.")

    return (focal_length_mm, working_distance_mm, fov_width_mm, fov_height_mm,
            sensor_w_mm_long, sensor_h_mm_short)


# ==============================================================================
# MAIN CLASS
# ==============================================================================

class DICSpeckleCalculator:
    def __init__(self,
                 resolution_width_px, resolution_height_px,
                 pixel_pitch_um,
                 pixels_per_speckle_range=(3, 5),
                 focal_length_mm=None,
                 working_distance_mm=None,
                 fov_width_mm=None):
        """
        Provide image resolution (px) and pixel pitch (µm).
        Provide EXACTLY TWO of:
          - focal_length_mm (f)
          - working_distance_mm (s)
          - fov_width_mm (FOV along the longest sensor side)

        FOV_height is automatically calculated from sensor aspect ratio.
        """
        self.resolution_width_px = int(resolution_width_px)
        self.resolution_height_px = int(resolution_height_px)
        self.pixel_pitch_um = float(pixel_pitch_um)
        self.pixels_per_speckle_range = tuple(pixels_per_speckle_range)

        # Eq. (1) Derive sensor size from resolution × pitch
        self.sensor_w_mm_native = (self.resolution_width_px * self.pixel_pitch_um) / 1000.0
        self.sensor_h_mm_native = (self.resolution_height_px * self.pixel_pitch_um) / 1000.0

        # Solve optics with automatic FOV_h calculation
        (self.focal_length_mm,
         self.working_distance_mm,
         self.fov_width_mm,
         self.fov_height_mm,
         self.sensor_w_mm_long,
         self.sensor_h_mm_short) = solve_optics_with_auto_fov(
            self.sensor_w_mm_native, self.sensor_h_mm_native,
            focal_length_mm=focal_length_mm,
            working_distance_mm=working_distance_mm,
            fov_width_mm=fov_width_mm
        )

        # Proceed with remaining calculations
        self.calculate()

    def calculate(self):
        """Calculate imaging scale and speckle sizes"""
        # Eq. (4) m = sensor_w / FOV_w
        self.magnification = self.sensor_w_mm_long / self.fov_width_mm

        # Pixel size on object and image scale (Eqs. (7) and (8))
        self.pixel_size_mm = calculate_pixel_size_on_object(self.pixel_pitch_um, self.magnification)
        self.image_scale_px_per_mm = calculate_image_scale(self.magnification, self.pixel_pitch_um)

        # Speckle sizes (Eq. (9))
        self.speckle_size_min_mm = calculate_speckle_size(self.pixels_per_speckle_range[0], self.image_scale_px_per_mm)
        self.speckle_size_max_mm = calculate_speckle_size(self.pixels_per_speckle_range[1], self.image_scale_px_per_mm)
        self.speckle_size_target_mm = 0.5 * (self.speckle_size_min_mm + self.speckle_size_max_mm)

        # Subset estimate
        speckles_per_subset = 5
        self.subset_size_min_px = speckles_per_subset * self.pixels_per_speckle_range[0]
        self.subset_size_max_px = speckles_per_subset * self.pixels_per_speckle_range[1]

        # Technique suggestion
        self.recommended_techniques = recommend_technique(self.speckle_size_target_mm)

    def _calculate_working_distances(self):
        """Calculate working distances using different optical models"""
        s_thin = self.focal_length_mm * (1.0 + 1.0 / self.magnification)
        s_thick = s_thin * 1.02
        s_bfd_65 = s_thin - (self.focal_length_mm * 0.65)
        s_bfd_70 = s_thin - (self.focal_length_mm * 0.70)
        s_bfd_75 = s_thin - (self.focal_length_mm * 0.75)

        return {
            'thin_lens': s_thin,
            'thick_lens': s_thick,
            'bfd_65': s_bfd_65,
            'bfd_70': s_bfd_70,
            'bfd_75': s_bfd_75,
        }

    def print_results(self):
        print("=" * 90)
        print("DIC SPECKLE PATTERN CALCULATOR - RESULTS")
        print("=" * 90)

        print("\n[INPUT]")
        print(f"  Resolution (px): {self.resolution_width_px} × {self.resolution_height_px}")
        print(f"  Pixel pitch:     {self.pixel_pitch_um:.3f} µm")
        print(f"  Pixels per speckle: {self.pixels_per_speckle_range[0]}-{self.pixels_per_speckle_range[1]} px")

        print("\n[DERIVED SENSOR]")
        print(f"  Sensor (native): {self.sensor_w_mm_native:.3f} × {self.sensor_h_mm_native:.3f} mm")
        print(f"  Sensor (long/short): {self.sensor_w_mm_long:.3f} × {self.sensor_h_mm_short:.3f} mm")
        print(f"  Aspect ratio: {self.sensor_w_mm_long / self.sensor_h_mm_short:.4f}")

        # Get working distances from different models
        wd = self._calculate_working_distances()

        print("\n[OPTICS]")
        print(f"  Focal length f:           {self.focal_length_mm:.3f} mm")
        print(f"  Field of view (W×H):      {self.fov_width_mm:.3f} × {self.fov_height_mm:.3f} mm")
        print(f"  Magnification m:          {self.magnification:.6f} ×")
        print(f"\n  Working distance s (optical models):")
        print(f"    Thin Lens (Principal Plane):     {wd['thin_lens']:.3f} mm")
        print(f"    Thick Lens (+2%):                {wd['thick_lens']:.3f} mm")
        print(f"    BFD Correction (65% of f):       {wd['bfd_65']:.3f} mm")
        print(f"    BFD Correction (70% of f):       {wd['bfd_70']:.3f} mm")
        print(f"    BFD Correction (75% of f):       {wd['bfd_75']:.3f} mm")

        print("\n[IMAGING SCALE]")
        print(f"  Pixel size on object: {self.pixel_size_mm * 1000:.2f} µm")
        print(f"  Image scale:          {self.image_scale_px_per_mm:.2f} px/mm")

        print("\n[SPECKLE PATTERN]")
        print(f"  Speckle size range:   {self.speckle_size_min_mm:.3f} – {self.speckle_size_max_mm:.3f} mm")
        print(
            f"  Optimal speckle size: {self.speckle_size_target_mm:.3f} mm ({self.speckle_size_target_mm * 1000:.1f} µm)")
        print(f"  Recommended subset:   {self.subset_size_min_px} – {self.subset_size_max_px} px")

        print("\n[RECOMMENDED PAINTING TECHNIQUES]")
        for i, t in enumerate(self.recommended_techniques, 1):
            print(f"  {i}. {t}")
        print("\n" + "=" * 90)


# ==============================================================================
# EXAMPLES
# ==============================================================================

if __name__ == "__main__":
    # Example A: Given f=50mm & FOV_w=25mm -> compute s
    print("\n>>> Example A: f=50mm & FOV_width=25mm")
    calc_A = DICSpeckleCalculator(
        resolution_width_px=2448, resolution_height_px=2048,
        pixel_pitch_um=3.45,
        focal_length_mm=50.0, fov_width_mm=25.0,
        pixels_per_speckle_range=(3, 5)
    )
    calc_A.print_results()

    # Example B: Given f=200mm & s=300mm -> compute FOV_w (and FOV_h auto)
    print("\n>>> Example B: f=200mm & s=300mm")
    calc_B = DICSpeckleCalculator(
        resolution_width_px=2448, resolution_height_px=2048,
        pixel_pitch_um=3.45,
        focal_length_mm=200.0, working_distance_mm=300.0,
        pixels_per_speckle_range=(3, 5)
    )
    calc_B.print_results()

    # Example C: Given s=400mm & FOV_w=45mm -> compute f
    print("\n>>> Example C: s=400mm & FOV_width=45mm")
    calc_C = DICSpeckleCalculator(
        resolution_width_px=2448, resolution_height_px=2048,
        pixel_pitch_um=3.45,
        working_distance_mm=400.0, fov_width_mm=45.0,
        pixels_per_speckle_range=(3, 5)
    )
    calc_C.print_results()

    print("\n✓ Workflow complete (FOV_height auto-calculated from sensor aspect ratio).")