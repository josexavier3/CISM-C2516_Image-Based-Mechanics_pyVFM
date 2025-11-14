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

import sys
from pathlib import Path
import numpy as np
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QComboBox, QLineEdit, QTextEdit, QGroupBox,
    QGridLayout, QSplitter, QFrame, QDoubleSpinBox
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QFontDatabase
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# ============================================================================
# GLOBAL FONT CONFIGURATION
# ============================================================================
FONT_FAMILY = "Inter"
FONT_FAMILY_MONO = "JetBrains Mono"
FONT_FALLBACK = "Segoe UI, Arial, sans-serif"
FONT_MONO_FALLBACK = "Consolas, Courier New, monospace"

FONT_SIZE_TITLE = 22
FONT_SIZE_SUBTITLE = 11
FONT_SIZE_HEADING = 12
FONT_SIZE_NORMAL = 10
FONT_SIZE_SMALL = 9
FONT_SIZE_MONO = 10

# ============================================================================
# GLOBAL COLOR CONFIGURATION
# ============================================================================
COLOR_ACCENT = "#F59E0B"
COLOR_ACCENT_HOVER = "#FBBF24"
COLOR_ACCENT_PRESSED = "#D97706"

COLOR_BUTTON_BG = "#52525B"
COLOR_BUTTON_HOVER = "#71717A"
COLOR_BUTTON_PRESSED = "#3F3F46"

COLOR_BG_MAIN = "#1a1a1a"
COLOR_BG_SECONDARY = "#242424"
COLOR_BG_INPUT = "#2d2d2d"
COLOR_BG_HEADER = "#1f1f1f"

COLOR_BORDER = "#404040"
COLOR_BORDER_FOCUS = "#F59E0B"
COLOR_SEPARATOR = "#404040"

COLOR_TEXT_PRIMARY = "#E5E7EB"
COLOR_TEXT_SECONDARY = "#9CA3AF"
COLOR_TEXT_MUTED = "#6B7280"

COLOR_VIZ_BLUR = "#EF4444"      # Red for blur line
COLOR_VIZ_LIMIT = "#10B981"     # Green for acceptable limit
COLOR_VIZ_SAFE = "#22D3EE"      # Cyan for safe zone
# ============================================================================


class BlurVisualizer(FigureCanvas):
    """Matplotlib canvas for motion blur visualization"""
    
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(6, 4), facecolor=COLOR_BG_SECONDARY)
        super().__init__(self.fig)
        self.setParent(parent)
        self.ax = self.fig.add_subplot(111, facecolor=COLOR_BG_SECONDARY)
        self.setup_plot()
        
    def setup_plot(self):
        """Initialize the plot"""
        self.ax.grid(True, alpha=0.3, color='white', linestyle='--')
        self.ax.set_xlabel('Exposure Time (ms)', color='white', fontsize=FONT_SIZE_NORMAL)
        self.ax.set_ylabel('Motion Blur (pixels)', color='white', fontsize=FONT_SIZE_NORMAL)
        self.ax.tick_params(colors='white', labelsize=FONT_SIZE_SMALL)
        
        for spine in self.ax.spines.values():
            spine.set_color('white')
            
        self.fig.tight_layout()
        
    def plot_blur_analysis(self, velocity, image_scale, max_blur, current_exposure, current_blur):
        """Plot motion blur vs exposure time"""
        self.ax.clear()
        self.setup_plot()
        
        # Generate exposure time range
        exposure_range = np.linspace(0.1, current_exposure * 3, 100)
        
        # Calculate blur for each exposure time
        blur_range = []
        for exp in exposure_range:
            velocity_mm_per_s = velocity / 60.0
            exposure_s = exp / 1000.0
            displacement = velocity_mm_per_s * exposure_s
            blur = displacement * image_scale
            blur_range.append(blur)
        
        # Plot acceptable limit (horizontal line)
        self.ax.axhline(y=max_blur, color=COLOR_VIZ_LIMIT, linestyle='--', 
                       linewidth=2, label=f'Acceptable Limit ({max_blur:.2f} px)')
        
        # Fill safe zone
        self.ax.fill_between(exposure_range, 0, max_blur, 
                            alpha=0.2, color=COLOR_VIZ_SAFE, label='Safe Zone')
        
        # Plot blur curve
        self.ax.plot(exposure_range, blur_range, color=COLOR_VIZ_BLUR, 
                    linewidth=2.5, label='Motion Blur')
        
        # Mark current operating point
        self.ax.plot(current_exposure, current_blur, 'o', 
                    color='#FCD34D', markersize=10, markeredgewidth=2,
                    markeredgecolor='white', label='Current Setting')
        
        # Add title
        self.ax.set_title('Motion Blur vs Exposure Time', color='white',
                         fontsize=FONT_SIZE_HEADING, fontweight='bold', pad=10)
        
        # Legend outside plot
        legend = self.ax.legend(
            loc='upper left',
            bbox_to_anchor=(1.02, 1),
            facecolor=COLOR_BG_SECONDARY,
            edgecolor='white',
            fontsize=FONT_SIZE_SMALL,
            labelcolor='white',
            framealpha=0.9
        )
        
        self.fig.tight_layout()
        self.draw()


class MotionBlurCalculatorGUI(QMainWindow):
    
    def __init__(self):
        super().__init__()
        
        self.sensor_dimensions = {
            "1/3.2\" (4.5mm)": 4.54,
            "1/2.5\" (5.8mm)": 5.76,
            "1/2.3\" (6.2mm)": 6.17,
            "2/3\" (8.8mm)": 8.8,
            "1.0\" (12.8mm)": 12.8,
            "1.1\" (13.2mm)": 13.2,
            "4/3\" (17.3mm)": 17.3,
        }
        
        self.sensor_resolutions = {
            4.54: 1024,
            5.76: 1024,
            6.17: 1024,
            8.8: 1280,
            12.8: 1920,
            13.2: 1920,
            17.3: 2048,
        }
        
        self.load_custom_fonts()
        self.init_ui()
        self.apply_stylesheet()
        
    def load_custom_fonts(self):
        """Load custom fonts"""
        font_dir = Path(__file__).parent / "fonts"
        if font_dir.exists():
            for font_file in font_dir.glob("*.ttf"):
                QFontDatabase.addApplicationFont(str(font_file))
            for font_file in font_dir.glob("*.otf"):
                QFontDatabase.addApplicationFont(str(font_file))
        
    def init_ui(self):
        """Initialize the user interface"""
        self.setWindowTitle("Motion Blur Calculator - DIC Tensile Testing")
        self.setGeometry(100, 100, 1200, 750)
        
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        
        # Header
        header = self.create_header()
        main_layout.addWidget(header)
        
        # Splitter for resizable panels
        splitter = QSplitter(Qt.Horizontal)
        
        # Left panel - Inputs and Visualization
        left_panel = self.create_left_panel()
        splitter.addWidget(left_panel)
        
        # Right panel - Results
        right_panel = self.create_right_panel()
        splitter.addWidget(right_panel)
        
        splitter.setSizes([600, 600])
        main_layout.addWidget(splitter)
        
    def create_header(self):
        """Create header section"""
        header = QFrame()
        header.setObjectName("header")
        header.setFixedHeight(80)
        
        layout = QVBoxLayout(header)
        layout.setContentsMargins(20, 10, 20, 10)
        
        title = QLabel("📸 Motion Blur Calculator")
        title.setObjectName("title")
        title.setFont(QFont(FONT_FAMILY, FONT_SIZE_TITLE, QFont.Bold))
        
        subtitle = QLabel("DIC Tensile Testing - Exposure Time Optimization")
        subtitle.setObjectName("subtitle")
        subtitle.setFont(QFont(FONT_FAMILY, FONT_SIZE_SUBTITLE))
        
        layout.addWidget(title)
        layout.addWidget(subtitle)
        
        return header
        
    def create_left_panel(self):
        """Create left panel with inputs and visualization"""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(15, 15, 15, 15)
        layout.setSpacing(15)
        
        # Input section
        input_group = self.create_input_section()
        layout.addWidget(input_group)
        
        # Calculate button
        self.calc_button = QPushButton("Calculate Motion Blur")
        self.calc_button.setObjectName("calcButton")
        self.calc_button.setFixedHeight(50)
        self.calc_button.setFont(QFont(FONT_FAMILY, FONT_SIZE_HEADING, QFont.Bold))
        self.calc_button.clicked.connect(self.calculate_blur)
        layout.addWidget(self.calc_button)
        
        # Visualization section
        viz_group = self.create_visualization_section()
        layout.addWidget(viz_group)
        
        return panel
        
    def create_input_section(self):
        """Create input fields section"""
        group = QGroupBox("Input Parameters")
        group.setFont(QFont(FONT_FAMILY, FONT_SIZE_HEADING, QFont.Bold))
        
        layout = QGridLayout()
        layout.setSpacing(12)
        layout.setContentsMargins(15, 20, 15, 15)
        
        row = 0
        
        # Sensor size
        sensor_label = QLabel("Sensor Size:")
        sensor_label.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(sensor_label, row, 0)
        
        self.sensor_combo = QComboBox()
        self.sensor_combo.addItems(list(self.sensor_dimensions.keys()))
        self.sensor_combo.setCurrentText("1.0\" (12.8mm)")
        self.sensor_combo.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(self.sensor_combo, row, 1)
        row += 1
        
        # Focal length
        focal_label = QLabel("Focal Length (mm):")
        focal_label.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(focal_label, row, 0)
        
        self.focal_input = QLineEdit("60")
        self.focal_input.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(self.focal_input, row, 1)
        row += 1
        
        # Working distance
        working_label = QLabel("Working Distance (mm):")
        working_label.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(working_label, row, 0)
        
        self.working_input = QLineEdit("500")
        self.working_input.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(self.working_input, row, 1)
        row += 1
        
        # Crosshead velocity
        velocity_label = QLabel("Crosshead Velocity (mm/min):")
        velocity_label.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(velocity_label, row, 0)
        
        self.velocity_input = QLineEdit("1.0")
        self.velocity_input.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(self.velocity_input, row, 1)
        row += 1
        
        # Exposure time
        exposure_label = QLabel("Exposure Time (ms):")
        exposure_label.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(exposure_label, row, 0)
        
        self.exposure_input = QLineEdit("10")
        self.exposure_input.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(self.exposure_input, row, 1)
        row += 1
        
        # Max acceptable blur
        blur_label = QLabel("Max Blur Criterion (px):")
        blur_label.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(blur_label, row, 0)
        
        self.blur_input = QLineEdit("0.01")
        self.blur_input.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(self.blur_input, row, 1)
        
        blur_hint = QLabel("Typically 0.01 px for DIC")
        blur_hint.setObjectName("hint")
        blur_hint.setFont(QFont(FONT_FAMILY, FONT_SIZE_SMALL))
        layout.addWidget(blur_hint, row + 1, 1)
        
        layout.setColumnStretch(1, 1)
        group.setLayout(layout)
        
        return group
        
    def create_visualization_section(self):
        """Create visualization display"""
        group = QGroupBox("Motion Blur Analysis")
        group.setFont(QFont(FONT_FAMILY, FONT_SIZE_HEADING, QFont.Bold))
        
        layout = QVBoxLayout()
        layout.setContentsMargins(10, 20, 10, 10)
        
        self.blur_canvas = BlurVisualizer()
        layout.addWidget(self.blur_canvas)
        
        group.setLayout(layout)
        return group
        
    def create_right_panel(self):
        """Create right panel with results"""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(15, 15, 15, 15)
        
        results_group = self.create_results_section()
        layout.addWidget(results_group)
        
        return panel
        
    def create_results_section(self):
        """Create results text display"""
        group = QGroupBox("Calculation Results")
        group.setFont(QFont(FONT_FAMILY, FONT_SIZE_HEADING, QFont.Bold))
        
        layout = QVBoxLayout()
        layout.setContentsMargins(15, 20, 15, 15)
        
        self.results_text = QTextEdit()
        self.results_text.setFont(QFont(FONT_FAMILY_MONO, FONT_SIZE_MONO))
        self.results_text.setReadOnly(True)
        
        welcome_msg = """
╔════════════════════════════════════════════════════════╗
║        Motion Blur Calculator for DIC                  ║
╚════════════════════════════════════════════════════════╝

Instructions:
  1. Select camera sensor size
  2. Enter lens focal length (mm)
  3. Enter working distance (mm)
  4. Set crosshead velocity (mm/min)
  5. Set exposure time (ms)
  6. Define max acceptable blur (typically 0.01 px)
  7. Click "Calculate Motion Blur"

Ready to analyze! 🚀
"""
        self.results_text.setText(welcome_msg)
        
        layout.addWidget(self.results_text)
        group.setLayout(layout)
        
        return group
        
    def calculate_blur(self):
        """Perform motion blur calculation"""
        try:
            # Get inputs
            sensor_str = self.sensor_combo.currentText()
            sensor_width_mm = self.sensor_dimensions[sensor_str]
            focal_length = float(self.focal_input.text())
            working_distance = float(self.working_input.text())
            velocity = float(self.velocity_input.text())
            exposure_time = float(self.exposure_input.text())
            max_blur = float(self.blur_input.text())
            
            # Validate
            if focal_length <= 0 or working_distance <= 0 or velocity < 0 or exposure_time <= 0:
                raise ValueError("Values must be positive")
            
            if working_distance <= focal_length:
                raise ValueError("Working distance must be greater than focal length")
            
            # Calculate FOV and image scale
            fov_width_mm = sensor_width_mm * ((working_distance - focal_length) / focal_length)
            
            closest_sensor = min(self.sensor_resolutions.keys(), 
                               key=lambda x: abs(x - sensor_width_mm))
            sensor_width_px = self.sensor_resolutions[closest_sensor]
            
            image_scale = sensor_width_px / fov_width_mm  # px/mm
            
            # Calculate motion blur
            velocity_mm_per_s = velocity / 60.0
            exposure_s = exposure_time / 1000.0
            displacement_mm = velocity_mm_per_s * exposure_s
            blur_px = displacement_mm * image_scale
            
            # Calculate max allowable exposure
            max_displacement_mm = max_blur / image_scale
            max_exposure_s = max_displacement_mm / velocity_mm_per_s
            max_exposure_ms = max_exposure_s * 1000.0
            
            # Check if acceptable
            is_acceptable = blur_px <= max_blur
            
            # Format results
            status_symbol = "✓" if is_acceptable else "✗"
            status_text = "ACCEPTABLE" if is_acceptable else "UNACCEPTABLE"
            
            if is_acceptable:
                margin = (max_blur - blur_px) / max_blur * 100
                status_detail = f"Safety margin: {margin:.1f}%"
            else:
                excess = (blur_px / max_blur - 1) * 100
                status_detail = f"Exceeds limit by {excess:.1f}%\n  Reduce exposure to ≤ {max_exposure_ms:.2f} ms"
            
            results = f"""
╔════════════════════════════════════════════════════════╗
║           MOTION BLUR ANALYSIS RESULTS                 ║
╚════════════════════════════════════════════════════════╝

📷 OPTICAL SETUP
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Sensor:             {sensor_str}
  Focal length:       {focal_length:.1f} mm
  Working distance:   {working_distance:.1f} mm
  Field of view:      {fov_width_mm:.2f} mm
  Sensor resolution:  {sensor_width_px} px
  Image scale:        {image_scale:.3f} px/mm
                      {1/image_scale:.4f} mm/px

🔬 TEST CONDITIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Crosshead velocity: {velocity:.3f} mm/min
                      {velocity/60:.4f} mm/s
  Exposure time:      {exposure_time:.2f} ms

📊 MOTION BLUR ASSESSMENT
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Calculated blur:    {blur_px:.5f} px
  Blur criterion:     ≤ {max_blur:.2f} px
  Max safe exposure:  {max_exposure_ms:.2f} ms

✓/✗ RESULT
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  {status_symbol} {status_text}
  {status_detail}

💡 RECOMMENDATIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  • Use exposure ≤ {max_exposure_ms:.2f} ms for this velocity
  • Consider increasing lighting if exposure too short
  • For higher velocities, reduce exposure proportionally
"""
            
            self.results_text.setText(results)
            
            # Update visualization
            self.blur_canvas.plot_blur_analysis(
                velocity, image_scale, max_blur, exposure_time, blur_px
            )
            
        except ValueError as e:
            error_msg = f"""
╔════════════════════════════════════════════════════════╗
║                   ⚠️  ERROR                            ║
╚════════════════════════════════════════════════════════╝

{str(e)}

Please check your inputs:
  • All values must be positive numbers
  • Working distance must be > focal length
  • Use decimal point (.) for decimals
"""
            self.results_text.setText(error_msg)
            
        except Exception as e:
            self.results_text.setText(f"Error: {str(e)}")
            
    def apply_stylesheet(self):
        """Apply dark theme stylesheet"""
        stylesheet = f"""
        QMainWindow {{
            background-color: {COLOR_BG_MAIN};
        }}
        QWidget {{
            background-color: {COLOR_BG_MAIN};
            color: {COLOR_TEXT_PRIMARY};
            font-family: {FONT_FAMILY}, {FONT_FALLBACK};
        }}
        QFrame#header {{
            background-color: {COLOR_BG_HEADER};
            border-bottom: 2px solid {COLOR_ACCENT};
        }}
        QLabel#title {{
            color: {COLOR_ACCENT};
        }}
        QLabel#subtitle {{
            color: {COLOR_TEXT_SECONDARY};
        }}
        QLabel#hint {{
            color: {COLOR_TEXT_MUTED};
            font-style: italic;
        }}
        QGroupBox {{
            border: 2px solid {COLOR_BORDER};
            border-radius: 8px;
            margin-top: 12px;
            padding-top: 10px;
            color: {COLOR_ACCENT};
            font-family: {FONT_FAMILY}, {FONT_FALLBACK};
        }}
        QGroupBox::title {{
            subcontrol-origin: margin;
            left: 15px;
            padding: 0 5px;
        }}
        QLineEdit, QComboBox {{
            background-color: {COLOR_BG_INPUT};
            border: 1px solid {COLOR_BORDER};
            border-radius: 5px;
            padding: 8px;
            color: {COLOR_TEXT_PRIMARY};
            selection-background-color: {COLOR_ACCENT};
            font-family: {FONT_FAMILY}, {FONT_FALLBACK};
        }}
        QLineEdit:focus, QComboBox:focus {{
            border: 2px solid {COLOR_BORDER_FOCUS};
        }}
        QComboBox::drop-down {{
            border: none;
            width: 30px;
        }}
        QComboBox::down-arrow {{
            image: none;
            border-left: 5px solid transparent;
            border-right: 5px solid transparent;
            border-top: 5px solid {COLOR_TEXT_PRIMARY};
            margin-right: 10px;
        }}
        QComboBox QAbstractItemView {{
            background-color: {COLOR_BG_INPUT};
            border: 1px solid {COLOR_BORDER};
            selection-background-color: {COLOR_ACCENT};
            color: {COLOR_TEXT_PRIMARY};
            font-family: {FONT_FAMILY}, {FONT_FALLBACK};
        }}
        QPushButton#calcButton {{
            background-color: {COLOR_BUTTON_BG};
            border: none;
            border-radius: 8px;
            color: {COLOR_TEXT_PRIMARY};
            padding: 10px;
            font-family: {FONT_FAMILY}, {FONT_FALLBACK};
        }}
        QPushButton#calcButton:hover {{
            background-color: {COLOR_BUTTON_HOVER};
        }}
        QPushButton#calcButton:pressed {{
            background-color: {COLOR_BUTTON_PRESSED};
        }}
        QTextEdit {{
            background-color: {COLOR_BG_INPUT};
            border: 1px solid {COLOR_BORDER};
            border-radius: 5px;
            padding: 10px;
            color: {COLOR_TEXT_PRIMARY};
            font-family: {FONT_FAMILY_MONO}, {FONT_MONO_FALLBACK};
        }}
        QSplitter::handle {{
            background-color: {COLOR_SEPARATOR};
            width: 2px;
        }}
        QSplitter::handle:hover {{
            background-color: {COLOR_ACCENT};
        }}
        """
        self.setStyleSheet(stylesheet)


def main():
    app = QApplication(sys.argv)
    window = MotionBlurCalculatorGUI()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
