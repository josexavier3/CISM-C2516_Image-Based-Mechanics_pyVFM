# Field of View (FOV) Calculator

A tool for calculating the field of view in Digital Image Correlation (DIC) experimental setups. Helps optimize camera and lens selection for your testing requirements.

## 🎯 Purpose

When setting up a DIC system, you need to ensure:
- Your specimen fits within the camera's field of view
- You have adequate resolution for strain measurements
- The magnification is appropriate for your application

This calculator eliminates guesswork by computing exact FOV dimensions based on your optical setup.

## ✨ Features

- **Calculate FOV dimensions** (width, height, diagonal) in mm
- **Determine optical magnification** 
- **Recommend Region of Interest (ROI)** for DIC analysis (90-95% of FOV)
- **Interactive visualization** showing FOV and recommended ROI boundaries
- **Support for common sensor sizes** (1/3.2", 1/2.5", 1/2.3", 2/3", 1", 1.1", 4/3")
- **Modern GUI interface** with dark theme

## 📸 Screenshots

### GUI Interface
The tool features a split-panel interface:
- **Left panel:** Input parameters and FOV visualization
- **Right panel:** Detailed calculation results

## 🚀 Quick Start

### GUI Version (Recommended)

```bash
python fov_calculator_gui_pyside6.py
```

### Command-Line Version

```bash
python fov_calculator.py
```

Edit the parameters at the bottom of `fov_calculator.py`:
```python
focal_length = 50        # mm
sensor_size = 1.0        # inches
standoff_distance = 500  # mm
```

## 📋 Requirements

```bash
pip install numpy matplotlib PySide6
```

## 📝 Input Parameters

| Parameter | Description | Typical Values |
|-----------|-------------|----------------|
| **Sensor Size** | Camera sensor diagonal size | 1/3.2", 1/2.5", 1/2.3", 2/3", 1", 1.1", 4/3" |
| **Focal Length** | Lens focal length | 16, 25, 35, 50, 75, 100 mm |
| **Stand-off Distance** | Distance from lens to specimen | 300-1000 mm |

## 📊 Output Information

The calculator provides:

### Optical Setup
- Sensor physical dimensions (mm)
- Field of view dimensions (width × height, mm)
- FOV diagonal (mm)
- Magnification factor

### DIC Recommendations
- Recommended ROI width range (90-95% of FOV)
- Recommended ROI height range (90-95% of FOV)

### Visualization
- Orange rectangle: Total field of view
- Green dashed rectangle: Recommended ROI (95% of FOV)
- Crosshairs: Optical center

## 🎨 Customization

The GUI supports customization through global configuration at the top of `fov_calculator_gui_pyside6.py`:

### Fonts
```python
FONT_FAMILY = "Inter"              # Main UI font
FONT_FAMILY_MONO = "JetBrains Mono"  # Results display
```

### Colors
```python
COLOR_ACCENT = "#F59E0B"           # Orange accent (headers, borders)
COLOR_BUTTON_BG = "#52525B"        # Grey button
COLOR_VIZ_FOV = "#F59E0B"          # FOV rectangle color
COLOR_VIZ_ROI = "#06B6D4"          # ROI rectangle color
```

## 💡 Usage Tips

### Choosing Stand-off Distance
- **Closer distance** → Smaller FOV, higher magnification, better resolution
- **Farther distance** → Larger FOV, lower magnification, more specimen coverage

### Optimal Setup
1. Start with your specimen dimensions
2. Add 10-15% margin for the FOV
3. Calculate required stand-off distance
4. Verify working distance allows adequate lighting

### Common Issues
- **FOV too small:** Increase stand-off distance or use shorter focal length
- **FOV too large:** Decrease stand-off distance or use longer focal length
- **Resolution insufficient:** Use larger sensor or move camera closer

## 📐 Theory

The FOV is calculated using thin lens equations:

```
FOV_width = Sensor_width × (Working_distance - Focal_length) / Focal_length
```

Where:
- `Sensor_width`: Physical sensor dimension (mm)
- `Working_distance`: Lens-to-specimen distance (mm)
- `Focal_length`: Lens focal length (mm)

Magnification:
```
M = Sensor_width / FOV_width
```

## 🔗 Related Tools

- **[Motion Blur Calculator](../motion_blur_calculator/)** - Optimize exposure time for your test velocity

## 📚 References

1. Sutton, M. A., Orteu, J. J., & Schreier, H. (2009). *Image correlation for shape, motion and deformation measurements: basic concepts, theory and applications*. Springer Science & Business Media.

2. Pan, B., Qian, K., Xie, H., & Asundi, A. (2009). Two-dimensional digital image correlation for in-plane displacement and strain measurement: a review. *Measurement science and technology*, 20(6), 062001.

## 🐛 Issues

Report issues to: jmc.xavier@fct.unl.pt

## 📄 License

GNU General Public License v3.0 (GPLv3)

---

**Part of the CISM C2516 Image-Based Mechanics Course Materials**  
José Xavier - Universidade NOVA de Lisboa
