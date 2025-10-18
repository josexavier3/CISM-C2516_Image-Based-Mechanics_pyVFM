# DIC Speckle Size Calculator

A comprehensive tool for calculating optimal speckle sizes in Digital Image Correlation (DIC) experimental setups. Helps optimize speckle pattern preparation and ensure adequate resolution for strain measurements.

## Purpose

When preparing a specimen for DIC, you need to ensure:
- Speckle size is appropriate for your camera resolution and optics
- Speckles are large enough to be resolved but small enough for measurement accuracy
- The painting technique matches your target speckle dimensions
- Subset sizes provide adequate texture for correlation

This calculator eliminates guesswork by computing exact speckle sizes based on your optical setup and camera parameters.

## Features

- **Calculate optimal speckle sizes** in mm based on optical parameters
- **Determine imaging scale** (px/mm) from magnification
- **Recommend painting techniques** for achieving target speckle sizes
- **Estimate subset dimensions** for DIC correlation windows
- **Support for common cameras** (Allied Vision, Basler, FLIR with factory defaults)
- **Flexible optical solver** - provide any 2 of 3 parameters:
  - Lens focal length
  - Working distance (object distance)
  - Field of view (width)
- **Multiple working distance models** (thin lens, thick lens, back focal distance corrections)
- **Modern GUI interface** with dark theme and live visualization

## Screenshots

### GUI Interface
The tool features a split-panel layout:
- **Left panel:** Camera parameters, optical parameters, and calculation button
- **Right panel:** Detailed calculation results with speckle recommendations

```
┌─────────────────────────────────────────────────────────────┐
│         🎯 DIC Speckle Pattern Calculator                  │
│         Optimal Speckle Size for DIC Measurements          │
├───────────────────────────┬──────────────────────────────────┤
│  Camera Parameters        │                                  │
│  ────────────────────────   Calculation Results             │
│  Camera: [Dropdown▼]      │  ┌──────────────────────────┐  │
│  Resolution: 4512 × 3211  │  │ 📷 CAMERA SETUP          │  │
│  Pixel Size: 2.74 µm      │  │ Sensor: 13.2 × 8.8 mm   │  │
│                            │  │ Resolution: 4512 × 3211 │  │
│  Optical Parameters        │  │ Pixel size: 2.74 µm     │  │
│  ────────────────────────   │                          │  │
│  ☑ Focal Length: 50 mm    │  │ 📏 OPTICS                │  │
│  ☑ Object Distance: 300   │  │ Focal length: 50 mm      │  │
│  ☐ Field of View: [──]    │  │ Distance: 300 mm         │  │
│                            │  │ FOV: 25.0 × 16.7 mm     │  │
│  [Calculate Speckle Size]  │  │ Magnification: 0.528×    │  │
│                            │  │                          │  │
│                            │  │ 🎯 SPECKLE SIZE          │  │
│                            │  │ Range: 0.057 – 0.095 mm │  │
│                            │  │ Target: 0.076 mm         │  │
│                            │  │                          │  │
│                            │  │ 🎨 TECHNIQUES            │  │
│                            │  │ • Spray paint (fine)     │  │
│                            │  └──────────────────────────┘  │
└───────────────────────────┴──────────────────────────────────┘
```

## Quick Start

### GUI Version (Recommended)

```bash
python -m dic_speckle.gui
```

Or if structured as a package:
```bash
pip install -e .
dic-calculator
```

### Command-Line Version

```bash
python src/dic_speckle/core.py
```

By default, this runs three example scenarios and prints results to console.

### Edit Parameters

For CLI use, edit the example section at the bottom of `src/dic_speckle/core.py`:

```python
calc = DICSpeckleCalculator(
    resolution_width_px=2448,
    resolution_height_px=2048,
    pixel_pitch_um=3.45,
    focal_length_mm=50.0,        # Provide two of three:
    fov_width_mm=25.0,         # - working_distance_mm (s)
    #working_distance_mm=300.0,   # - focal_length_mm (f)
    pixels_per_speckle_range=(3, 5)
)
calc.print_results()
```

## Requirements

### Core Dependencies
```bash
pip install numpy
```

### GUI Dependencies (Optional)
```bash
pip install PySide6
```

### Development
```bash
pip install pytest black pylint
```

## Input Parameters

| Parameter | Description | Typical Values | Units |
|-----------|-------------|----------------|-------|
| **Resolution (H × V)** | Image resolution in pixels | 2048–4512 | px |
| **Pixel Pitch** | Sensor pixel size | 2.74–5.4 | µm |
| **Focal Length (f)** | Lens focal length | 16–200 | mm |
| **Working Distance (s)** | Object plane to lens principal plane | 300–1000 | mm |
| **Field of View Width** | Coverage on specimen (calculated) | 10–500 | mm |
| **Pixels per Speckle** | Desired speckle resolution | 3–7 | px |

## Output Information

The calculator provides comprehensive results:

### Optical Setup
- Sensor physical dimensions (mm)
- Field of view dimensions (width × height, mm)
- Magnification factor (m)
- Working distances under different optical models:
  - Thin lens (principal plane reference)
  - Thick lens (±2% correction)
  - Back focal distance (BFD) corrections (65%, 70%, 75% of focal length)

### Imaging Scale
- Pixel size on object (µm)
- Image scale (px/mm)

### Speckle Pattern Recommendations
- **Minimum speckle size** (lower pixels_per_speckle bound)
- **Maximum speckle size** (upper pixels_per_speckle bound)
- **Target speckle size** (midpoint)
- **Recommended painting techniques** that match target size:
  - Fine airbrush (0.3mm nozzle): 0.05–0.30 mm
  - Standard airbrush (0.5mm): 0.20–0.80 mm
  - Spray paint (fine): 0.40–1.50 mm
  - Spray paint (coarse): 0.80–3.00 mm
  - Speckle stamp: 0.30–1.50 mm
  - Toner powder: 0.01–0.10 mm
  - Lithography/etching: 0.001–0.050 mm

### DIC Recommendations
- **Subset size range** (speckles_per_subset × pixels_per_speckle)
- For typical DIC: 5 speckles per subset → 15–25 px subsets

## Customization

### Adding Custom Cameras

```python
CAMERA_DATABASE = {
    "My Custom Camera": {
        "sensor_width_mm": 17.3,
        "sensor_height_mm": 13.0,
        "pixel_size_um": 4.5,
        "resolution_h": 3840,
        "resolution_v": 2880,
        "sensor": "Sony IMX-XXXX",
    },
    # ... existing cameras ...
}
```

### Modifying Painting Techniques

```python
PAINTING_TECHNIQUES = {
    'My technique': (min_size_mm, max_size_mm),
    'Fine airbrush (0.3mm nozzle)': (0.05, 0.30),
    # ...
}
```

### GUI Theme Customization

```python
COLOR_ACCENT = "#F59E0B"           # Orange (headers, highlights)
COLOR_BG_MAIN = "#1a1a1a"          # Dark background
COLOR_TEXT_PRIMARY = "#E5E7EB"     # Light text
```

## Usage Tips

### Choosing Speckle Size

**Too small** (< 3 px):
- Risk: Speckles not resolved by camera
- Result: Loss of correlation signal

**Optimal** (3–7 px):
- Best correlation quality
- Typical range: 3–5 px

**Too large** (> 7 px):
- Risk: Reduced spatial resolution
- Result: Loss of strain measurement detail

### Workflow

1. **Select your camera** or input custom specifications
2. **Provide two optical parameters:**
   - If you know focal length and working distance → calculates FOV
   - If you know focal length and desired FOV → calculates working distance
   - If you know working distance and desired FOV → calculates focal length
3. **Set pixels per speckle range** (typically 3–5)
4. **Review recommended speckle size**
5. **Select painting technique** from the recommendations
6. **Apply pattern** to your specimen

### Common Scenarios

#### Scenario 1: Fixed Lens, Variable Distance
- Known: f = 50mm, want 0.1 mm speckles
- Solution: Adjust working distance until target size achieved
- Use case: Specimen orientation flexibility

#### Scenario 2: Fixed Position, Variable Magnification
- Known: s = 500mm, want to cover 50×50 mm specimen
- Solution: Adjust focal length or use different lens
- Use case: Fixed support structure

#### Scenario 3: Fixed Camera, New Application
- Known: All camera parameters, resolution fixed
- Solution: Choose focal length and working distance for application
- Use case: New test setup

### Troubleshooting

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| Speckle size too small | High magnification or small pixels | Increase working distance or use shorter focal length |
| Speckle size too large | Low magnification | Decrease working distance or use longer focal length |
| FOV too small | Magnification too high | Increase working distance |
| FOV too large | Magnification too low | Decrease working distance |
| Pattern application difficult | No matching technique | Consider custom approach or adjust target size |

## Theory

### Magnification

From thin lens equation with object at distance s and focal length f:

```
m = f / (s - f)
```

Where:
- m = magnification (image/object ratio)
- f = focal length (mm)
- s = working distance (mm)

Equivalently, using field of view:

```
m = Sensor_width / FOV_width
```

### Pixel Size on Object

The size of each camera pixel when projected onto the object:

```
p_obj = pixel_pitch / (m × 1000)
```

Where:
- p_obj = pixel size on object (mm)
- pixel_pitch = sensor pixel size (µm)
- m = magnification (dimensionless)

### Image Scale

Converts pixel coordinates to object coordinates:

```
scale = m × 1000 / pixel_pitch    [px/mm]
```

### Speckle Size in Pixels vs. Millimeters

If a speckle occupies N pixels on the image:

```
speckle_size [mm] = N [px] / scale [px/mm]
```

Typically, N = 3–5 for optimal DIC correlation.

## Project Structure

```
dic_speckle_calculator/
├── README.md
├── speckle_size_calculator.py
|── speckle_size_calculator_gui

```
## Issues & Support

Report issues to: jmc.xavier@fct.unl.pt

For bug reports, please include:
- Operating system and Python version
- Input parameters that caused the issue
- Full error message and traceback

## License

GNU General Public License v3.0 (GPLv3)

---

**Part of the CISM C2516 Course Materials**  
*"Image-Based Mechanics: An Overview of Experimental and Numerical Approaches"*  

Developed by: **José Xavier**  
Universidade NOVA de Lisboa, NOVA FCT, UNIDEMI  
https://userweb.fct.unl.pt/~jmc.xavier/index.html

Course information: https://cism.it/en/activities/courses/C2516/