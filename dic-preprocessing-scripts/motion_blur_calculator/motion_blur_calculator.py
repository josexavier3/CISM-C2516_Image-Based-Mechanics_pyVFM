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
"""

import numpy as np

def calculate_image_scale(sensor_width_mm, focal_length_mm, working_distance_mm):
    """
    Calculate image scale (pixels per mm on specimen)
    
    Parameters:
    -----------
    sensor_width_mm : float
        Physical width of sensor in mm
    focal_length_mm : float
        Lens focal length in mm
    working_distance_mm : float
        Distance from lens to specimen in mm
    
    Returns:
    --------
    float : Image scale in px/mm
    """
    # Field of view on specimen (mm)
    fov_width_mm = sensor_width_mm * (working_distance_mm - focal_length_mm) / focal_length_mm
    
    # Assuming standard sensor resolutions based on sensor size
    # These are typical resolutions - adjust if needed
    sensor_resolutions = {
        12.8: 1920,   # 1" sensor → ~2MP
        13.2: 1920,   # 1.1" sensor → ~2MP
        8.8: 1280,    # 2/3" sensor → ~1.3MP
        6.17: 1024,   # 1/2.3" sensor → ~1MP
    }
    
    # Find closest sensor size
    closest_sensor = min(sensor_resolutions.keys(), 
                        key=lambda x: abs(x - sensor_width_mm))
    sensor_width_px = sensor_resolutions[closest_sensor]
    
    # Image scale: pixels per mm
    image_scale = sensor_width_px / fov_width_mm
    
    return image_scale, fov_width_mm, sensor_width_px


def calculate_motion_blur(velocity_mm_per_min, image_scale_px_per_mm, 
                          exposure_time_ms):
    """
    Calculate motion blur during exposure
    
    Parameters:
    -----------
    velocity_mm_per_min : float
        Specimen velocity in mm/min
    image_scale_px_per_mm : float
        Image scale in pixels per mm
    exposure_time_ms : float
        Exposure time in milliseconds
    
    Returns:
    --------
    float : Motion blur in pixels
    """
    # Convert velocity to mm/s
    velocity_mm_per_s = velocity_mm_per_min / 60.0
    
    # Convert exposure time to seconds
    exposure_time_s = exposure_time_ms / 1000.0
    
    # Calculate displacement during exposure
    displacement_mm = velocity_mm_per_s * exposure_time_s
    
    # Convert to pixels
    blur_px = displacement_mm * image_scale_px_per_mm
    
    return blur_px


def calculate_max_exposure(velocity_mm_per_min, image_scale_px_per_mm, 
                           max_blur_px=0.01):
    """
    Calculate maximum allowable exposure time
    
    Parameters:
    -----------
    velocity_mm_per_min : float
        Specimen velocity in mm/min
    image_scale_px_per_mm : float
        Image scale in pixels per mm
    max_blur_px : float
        Maximum acceptable blur in pixels (default 0.01)
    
    Returns:
    --------
    float : Maximum exposure time in milliseconds
    """
    # Convert velocity to mm/s
    velocity_mm_per_s = velocity_mm_per_min / 60.0
    
    # Maximum displacement in mm
    max_displacement_mm = max_blur_px / image_scale_px_per_mm
    
    # Maximum exposure time in seconds
    max_exposure_s = max_displacement_mm / velocity_mm_per_s
    
    # Convert to milliseconds
    max_exposure_ms = max_exposure_s * 1000.0
    
    return max_exposure_ms


def motion_blur_analysis(velocity_mm_per_min, sensor_width_mm, 
                        focal_length_mm, working_distance_mm,
                        exposure_time_ms, max_blur_px=0.01):
    """
    Complete motion blur analysis
    """
    
    # Calculate image scale
    image_scale, fov_width, sensor_width_px = calculate_image_scale(
        sensor_width_mm, focal_length_mm, working_distance_mm)
    
    # Calculate motion blur
    blur_px = calculate_motion_blur(velocity_mm_per_min, image_scale, 
                                    exposure_time_ms)
    
    # Calculate maximum allowable exposure
    max_exp_ms = calculate_max_exposure(velocity_mm_per_min, image_scale, 
                                        max_blur_px)
    
    # Check if acceptable
    is_acceptable = blur_px <= max_blur_px
    
    return {
        'image_scale_px_per_mm': image_scale,
        'fov_width_mm': fov_width,
        'sensor_width_px': sensor_width_px,
        'motion_blur_px': blur_px,
        'max_exposure_ms': max_exp_ms,
        'is_acceptable': is_acceptable,
        'blur_criterion_px': max_blur_px
    }


def print_analysis_report(result, velocity, exposure):
    """Print formatted analysis report"""
    print("=" * 70)
    print("MOTION BLUR ANALYSIS FOR DIC")
    print("=" * 70)
    
    print(f"\nOptical Setup:")
    print(f"  Field of View: {result['fov_width_mm']:.2f} mm")
    print(f"  Sensor Resolution: {result['sensor_width_px']} px")
    print(f"  Image Scale: {result['image_scale_px_per_mm']:.3f} px/mm")
    print(f"                = {1/result['image_scale_px_per_mm']:.4f} mm/px")
    
    print(f"\nTest Conditions:")
    print(f"  Crosshead Velocity: {velocity:.3f} mm/min = {velocity/60:.4f} mm/s")
    print(f"  Exposure Time: {exposure:.2f} ms")
    
    print(f"\nMotion Blur Assessment:")
    print(f"  Calculated Blur: {result['motion_blur_px']:.5f} px")
    print(f"  Blur Criterion: ≤ {result['blur_criterion_px']:.2f} px")
    print(f"  Maximum Safe Exposure: {result['max_exposure_ms']:.2f} ms")
    
    print(f"\nResult:")
    if result['is_acceptable']:
        print(f"  ✓ ACCEPTABLE - Blur within limits")
        margin = (result['blur_criterion_px'] - result['motion_blur_px']) / result['blur_criterion_px'] * 100
        print(f"  Safety margin: {margin:.1f}%")
    else:
        print(f"  ✗ UNACCEPTABLE - Reduce exposure time!")
        excess = (result['motion_blur_px'] / result['blur_criterion_px'] - 1) * 100
        print(f"  Exceeds limit by: {excess:.1f}%")
        print(f"  Recommended: Reduce exposure to ≤ {result['max_exposure_ms']:.2f} ms")
    
    print("=" * 70)


# Standard sensor sizes (diagonal in inches → width in mm)
SENSOR_SIZES = {
    '1/3.2"': 4.54,
    '1/2.5"': 5.76,
    '1/2.3"': 6.17,
    '2/3"': 8.8,
    '1"': 12.8,
    '1.1"': 13.2,
    '4/3"': 17.3,
}


if __name__ == "__main__":
    print("\n" + "="*70)
    print("MOTION BLUR CALCULATOR FOR TENSILE TESTING DIC")
    print("="*70)
    
    # ========== USER ADJUSTABLE PARAMETERS ==========
    
    # Test parameters
    velocity_mm_per_min = 1.0  # Crosshead velocity [mm/min]
    
    # Camera and lens parameters
    sensor_size = '1"'  # Choose from SENSOR_SIZES keys
    sensor_width_mm = SENSOR_SIZES[sensor_size]
    focal_length_mm = 60  # Lens focal length [mm]
    working_distance_mm = 500  # Distance to specimen [mm]
    
    # Exposure setting
    exposure_time_ms = 10  # Exposure time [ms]
    
    # Blur criterion
    max_blur_px = 0.01  # Maximum acceptable blur [px]
    
    # ================================================
    
    print(f"\nInput Parameters:")
    print(f"  Sensor: {sensor_size} ({sensor_width_mm} mm width)")
    print(f"  Focal Length: {focal_length_mm} mm")
    print(f"  Working Distance: {working_distance_mm} mm")
    print(f"  Crosshead Velocity: {velocity_mm_per_min} mm/min")
    print(f"  Exposure Time: {exposure_time_ms} ms")
    print()
    
    # Run analysis
    result = motion_blur_analysis(
        velocity_mm_per_min=velocity_mm_per_min,
        sensor_width_mm=sensor_width_mm,
        focal_length_mm=focal_length_mm,
        working_distance_mm=working_distance_mm,
        exposure_time_ms=exposure_time_ms,
        max_blur_px=max_blur_px
    )
    
    # Print results
    print_analysis_report(result, velocity_mm_per_min, exposure_time_ms)
    
    # Additional scenarios
    print("\n" + "="*70)
    print("TESTING DIFFERENT VELOCITIES")
    print("="*70)
    
    test_velocities = [0.1, 0.5, 1.0, 5.0, 10.0]  # mm/min
    
    print(f"\nFor exposure time = {exposure_time_ms} ms:")
    print(f"{'Velocity [mm/min]':>18} | {'Blur [px]':>12} | {'Max Exp [ms]':>14} | {'Status':>10}")
    print("-" * 70)
    
    for vel in test_velocities:
        res = motion_blur_analysis(vel, sensor_width_mm, focal_length_mm,
                                   working_distance_mm, exposure_time_ms, max_blur_px)
        status = "OK" if res['is_acceptable'] else "EXCEED"
        print(f"{vel:>18.2f} | {res['motion_blur_px']:>12.5f} | {res['max_exposure_ms']:>14.2f} | {status:>10}")
