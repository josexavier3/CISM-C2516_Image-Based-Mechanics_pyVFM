# Motion Blur Calculator

A tool for analyzing and preventing motion blur in Digital Image Correlation (DIC) measurements during tensile testing. Ensures optimal camera exposure settings for accurate strain measurements.

## 🎯 Purpose

During mechanical testing, specimen motion combined with camera exposure time causes motion blur, which can:
- Degrade DIC correlation quality
- Introduce measurement errors
- Reduce strain accuracy
- Cause correlation failures

This calculator determines the maximum allowable exposure time to keep motion blur below acceptable limits (typically ≤0.01 pixels for DIC).

## ✨ Features

- **Calculate motion blur** based on velocity and exposure time
- **Determine maximum safe exposure** for given test conditions
- **Assess acceptability** against DIC blur criteria
- **Interactive visualization** showing blur vs. exposure relationship
- **Safety margins** and recommendations
- **Multi-velocity analysis** for test planning
- **Modern GUI interface** with real-time calculations

## 📸 Screenshots

### GUI Interface
The tool features a split-panel interface:
- **Left panel:** Input parameters and blur analysis visualization
- **Right panel:** Detailed calculation results with pass/fail assessment

## 🚀 Quick Start

### GUI Version (Recommended)

```bash
python motion_blur_calculator_gui.py
```

### Command-Line Version

```bash
python motion_blur_calculator.py
```

Edit the parameters in `motion_blur_calculator.py`:
```python
velocity_mm_per_min = 1.0      # Crosshead velocity
sensor_size = '1"'              # Camera sensor
focal_length_mm = 60            # Lens focal length
working_distance_mm = 500       # Distance to specimen
exposure_time_ms = 10           # Exposure time
max_blur_px = 0.01              # Acceptable blur limit
```

## 📋 Requirements

```bash
pip install numpy matplotlib PySide6
```

## 📝 Input Parameters

| Parameter | Description | Typical Values | Units |
|-----------|-------------|----------------|-------|
| **Sensor Size** | Camera sensor diagonal | 1/3.2", 1/2.5", 2/3", 1", 1.1", 4/3" | inches |
| **Focal Length** | Lens focal length | 16, 25, 35, 50, 60, 75, 100 | mm |
| **Working Distance** | Lens to specimen distance | 300-1000 | mm |
| **Crosshead Velocity** | Test machine speed | 0.1-10 | mm/min |
| **Exposure Time** | Camera exposure duration | 1-50 | ms |
| **Max Blur Criterion** | Acceptable blur limit | 0.01 (DIC standard) | pixels |

## 📊 Output Information

The calculator provides:

### Optical Setup
- Field of view dimensions
- Sensor resolution
- Image scale (px/mm and mm/px)

### Test Conditions
- Crosshead velocity (mm/min and mm/s)
- Exposure time (ms)

### Motion Blur Assessment
- Calculated blur (pixels)
- Blur criterion threshold
- Maximum safe exposure time
- **Pass/Fail status** with safety margin or excess percentage

### Recommendations
- Suggested maximum exposure time
- Lighting considerations
- Velocity-dependent adjustments

### Visualization
- **Red curve:** Motion blur vs. exposure time
- **Green line:** Acceptable blur limit
- **Cyan zone:** Safe operating region
- **Yellow dot:** Current exposure setting

## 🎨 Customization

Customize the GUI by editing global configuration in `motion_blur_calculator_gui.py`:

### Fonts
```python
FONT_FAMILY = "Inter"
FONT_FAMILY_MONO = "JetBrains Mono"
```

### Colors
```python
COLOR_ACCENT = "#F59E0B"           # Orange accent
COLOR_BUTTON_BG = "#52525B"        # Grey button
COLOR_VIZ_BLUR = "#EF4444"         # Red - blur curve
COLOR_VIZ_LIMIT = "#10B981"        # Green - acceptable limit
COLOR_VIZ_SAFE = "#22D3EE"         # Cyan - safe zone
```

## 💡 Usage Tips

### Understanding the Results

#### ✓ ACCEPTABLE
```
Motion blur: 0.0087 px
Criterion: ≤ 0.01 px
Safety margin: 13.2%
```
Your settings are safe. You have 13% margin before reaching the limit.

#### ✗ UNACCEPTABLE
```
Motion blur: 0.0134 px
Criterion: ≤ 0.01 px
Exceeds limit by: 34.2%
Recommended: Reduce exposure to ≤ 7.46 ms
```
Reduce your exposure time to the recommended value or lower.

### Optimizing Your Setup

1. **Start with desired velocity** (based on material/standard)
2. **Calculate maximum exposure** using this tool
3. **Adjust lighting** if exposure is too short
4. **Verify correlation quality** in practice

### Common Scenarios

| Test Speed | Typical Max Exposure | Notes |
|------------|---------------------|-------|
| 0.1 mm/min | ~100 ms | Very slow tests, lighting flexible |
| 1.0 mm/min | ~10 ms | Standard rate, moderate lighting needed |
| 5.0 mm/min | ~2 ms | Fast tests, bright lighting required |
| 10 mm/min | ~1 ms | Very fast, high-power lighting essential |

### Troubleshooting

**Problem:** Calculated exposure too short for adequate lighting
**Solutions:**
- Increase lighting intensity
- Use brighter LED arrays
- Reduce test velocity (if possible)
- Use larger aperture lens
- Increase ISO (reduces SNR)

**Problem:** Camera can't achieve required exposure time
**Solutions:**
- Use faster camera with shorter minimum exposure
- Reduce test velocity
- Accept slightly higher blur (test correlation quality)

## 📐 Theory

### Motion Blur Calculation

```python
# Specimen displacement during exposure
velocity_mm_s = velocity_mm_min / 60
exposure_s = exposure_ms / 1000
displacement_mm = velocity_mm_s × exposure_s

# Convert to pixels using image scale
blur_px = displacement_mm × image_scale_px_per_mm
```

### Maximum Exposure Time

```python
max_displacement_mm = max_blur_px / image_scale
max_exposure_s = max_displacement_mm / velocity_mm_s
max_exposure_ms = max_exposure_s × 1000
```

### Why 0.01 Pixels?

The 0.01 pixel criterion comes from DIC accuracy requirements:
- DIC can measure displacements to ~0.01-0.05 px accuracy
- Motion blur should be ≤ measurement precision
- 0.01 px ensures blur doesn't dominate measurement uncertainty

## 🔬 Practical Example

### Setup
- Material: Aluminum alloy
- Test: Tensile at 1 mm/min
- Camera: 1" sensor, 1920×1080 px
- Lens: 50mm focal length
- Working distance: 500mm

### Results
- Image scale: 16.16 px/mm
- Maximum exposure: 9.92 ms
- **Recommendation:** Use 8-9 ms exposure

### Implementation
- Set camera exposure to 8 ms
- Verify adequate brightness
- Start test and confirm correlation quality

## 🔗 Related Tools

- **[FOV Calculator](../fov_calculator/)** - Plan your optical setup and field of view

## 📚 References

1. Sutton, M. A., et al. (2009). *Image correlation for shape, motion and deformation measurements*. Springer.

2. Pan, B. (2018). Digital image correlation for surface deformation measurement: historical developments, recent advances and future goals. *Measurement Science and Technology*, 29(8), 082001.

3. Lecompte, D., et al. (2006). Quality assessment of speckle patterns for digital image correlation. *Optics and Lasers in Engineering*, 44(11), 1132-1145.

## 🐛 Issues

Report issues to: jmc.xavier@fct.unl.pt

## 📄 License

GNU General Public License v3.0 (GPLv3)

---

**Part of the CISM C2516 Image-Based Mechanics Course Materials**  
José Xavier - Universidade NOVA de Lisboa
