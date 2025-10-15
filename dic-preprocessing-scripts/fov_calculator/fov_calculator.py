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

def calculate_fov(sensor_size_inches, focal_length_mm, standoff_distance_mm):
    """
    Calculate Field of View (FOV) for DIC setup
    
    Parameters:
    -----------
    sensor_size_inches : float
        Sensor diagonal size in inches (e.g., 1.1 for 1.1" sensor)
    focal_length_mm : float
        Lens focal length in millimeters
    standoff_distance_mm : float
        Distance from lens to object in millimeters
    
    Returns:
    --------
    dict : FOV dimensions in mm
    """
    
    # Standard sensor dimensions for common formats
    sensor_dimensions = {
        1/3.2: (4.54, 3.42),   # 1/3.2"
        1/2.5: (5.76, 4.29),   # 1/2.5"
        1/2.3: (6.17, 4.55),   # 1/2.3"
        2/3: (8.8, 6.6),       # 2/3"
        1.0: (12.8, 9.6),      # 1"
        1.1: (13.2, 8.8),      # 1.1"
        4/3: (17.3, 13.0),     # 4/3"
    }
    
    if sensor_size_inches in sensor_dimensions:
        L_CS_width, L_CS_height = sensor_dimensions[sensor_size_inches]
    else:
        print(f"Warning: Using approximate conversion for {sensor_size_inches}\" sensor")
        L_CS_width = sensor_size_inches * 25.4 * 0.6
        L_CS_height = L_CS_width * 0.75
    
    L_FL = focal_length_mm
    L_SOD = standoff_distance_mm
    
    # Calculate FOV using the formula
    L_FOV_width = L_CS_width * ((L_SOD - L_FL) / L_FL)
    L_FOV_height = L_CS_height * ((L_SOD - L_FL) / L_FL)
    
    # Calculate magnification
    magnification = L_CS_width / L_FOV_width
    
    # Calculate diagonal FOV
    L_FOV_diagonal = np.sqrt(L_FOV_width**2 + L_FOV_height**2)
    
    return {
        'sensor_width_mm': L_CS_width,
        'sensor_height_mm': L_CS_height,
        'fov_width_mm': L_FOV_width,
        'fov_height_mm': L_FOV_height,
        'fov_diagonal_mm': L_FOV_diagonal,
        'magnification': magnification,
        'focal_length_mm': L_FL,
        'standoff_distance_mm': L_SOD
    }


def print_results(result, sensor_size, focal_length, standoff_distance):
    """Print formatted results"""
    print("=" * 60)
    print("FOV CALCULATION FOR DIC SETUP")
    print("=" * 60)
    print(f"\nInput Parameters:")
    print(f"  Sensor size: {sensor_size}\"")
    print(f"  Sensor dimensions: {result['sensor_width_mm']:.2f} × {result['sensor_height_mm']:.2f} mm")
    print(f"  Focal length: {focal_length} mm")
    print(f"  Stand-off distance: {standoff_distance} mm")
    print(f"\nCalculated Field of View:")
    print(f"  FOV width: {result['fov_width_mm']:.2f} mm")
    print(f"  FOV height: {result['fov_height_mm']:.2f} mm")
    print(f"  FOV diagonal: {result['fov_diagonal_mm']:.2f} mm")
    print(f"  Magnification: {result['magnification']:.4f}×")
    print(f"\nFor DIC: ROI should be ~90-95% of FOV")
    print(f"  Recommended ROI width: {result['fov_width_mm'] * 0.9:.2f} - {result['fov_width_mm'] * 0.95:.2f} mm")
    print(f"  Recommended ROI height: {result['fov_height_mm'] * 0.9:.2f} - {result['fov_height_mm'] * 0.95:.2f} mm")
    print("=" * 60)


def interactive_fov_calculator():
    """Interactive FOV calculator"""
    print("\n" + "="*60)
    print("INTERACTIVE FOV CALCULATOR FOR DIC")
    print("="*60)
    
    sensor_size = float(input("\nEnter sensor size in inches (e.g., 1.1): "))
    focal_length = float(input("Enter focal length in mm (e.g., 50): "))
    standoff_distance = float(input("Enter stand-off distance in mm (e.g., 500): "))
    
    result = calculate_fov(sensor_size, focal_length, standoff_distance)
    print_results(result, sensor_size, focal_length, standoff_distance)


if __name__ == "__main__":
    # Example calculation
    focal_length = 50  # mm
    sensor_size = 1.  # inches (1.1" sensor); 1/3.2", 1/2.5", 1/2.3", 2/3", 1", 1.1", 4/3"
    standoff_distance = 500  # mm (example: 50 cm from object)
    
    result = calculate_fov(sensor_size, focal_length, standoff_distance)
    print_results(result, sensor_size, focal_length, standoff_distance)
    
    # Uncomment the line below for interactive mode
    # interactive_fov_calculator()
