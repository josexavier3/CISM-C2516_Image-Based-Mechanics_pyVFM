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
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QComboBox, QLineEdit, QTextEdit, QGroupBox,
    QGridLayout, QFrame, QCheckBox
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QFontDatabase

# Import the calculator module
from speckle_size_calculator import (
    DICSpeckleCalculator,
    PAINTING_TECHNIQUES,
    recommend_technique
)

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

# ============================================================================
# CAMERA DATABASE
# ============================================================================
CAMERA_DATABASE = {
    "Custom": {
        "sensor_width_mm": 12.8,
        "sensor_height_mm": 9.6,
        "pixel_size_um": 5.23,
        "resolution_h": 2448,
        "resolution_v": 2048,
    },
    "Allied Vision Alvium 1800 U-2050 (IMX541)": {
        "sensor_width_mm": 13.2,
        "sensor_height_mm": 8.8,
        "pixel_size_um": 2.74,
        "resolution_h": 4512,
        "resolution_v": 3211,
        "sensor": "Sony IMX541",
    },
    "Basler acA2440-75um (IMX392)": {
        "sensor_width_mm": 13.2,
        "sensor_height_mm": 8.8,
        "pixel_size_um": 5.4,
        "resolution_h": 2448,
        "resolution_v": 2048,
        "sensor": "Sony IMX392",
    },
    "FLIR Blackfly S BFS-U3-51S5M": {
        "sensor_width_mm": 8.8,
        "sensor_height_mm": 6.6,
        "pixel_size_um": 3.45,
        "resolution_h": 2448,
        "resolution_v": 2048,
        "sensor": "Sony IMX250",
    },
}


class DICSpeckleCalculatorGUI(QMainWindow):

    def __init__(self):
        super().__init__()
        self.calculator = None
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
        """Initialize UI"""
        self.setWindowTitle("DIC Speckle Pattern Calculator")
        self.setGeometry(100, 100, 1200, 700)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Header
        header = self.create_header()
        main_layout.addWidget(header)

        # Content area
        content = QWidget()
        content_layout = QHBoxLayout(content)
        content_layout.setContentsMargins(15, 15, 15, 15)
        content_layout.setSpacing(15)

        # Left panel - Inputs
        left_panel = self.create_left_panel()
        content_layout.addWidget(left_panel, 2)

        # Right panel - Results
        right_panel = self.create_right_panel()
        content_layout.addWidget(right_panel, 3)

        main_layout.addWidget(content)

    def create_header(self):
        """Create header"""
        header = QFrame()
        header.setObjectName("header")
        header.setFixedHeight(80)

        layout = QVBoxLayout(header)
        layout.setContentsMargins(20, 10, 20, 10)

        title = QLabel("🎯 DIC Speckle Pattern Calculator")
        title.setObjectName("title")
        title.setFont(QFont(FONT_FAMILY, FONT_SIZE_TITLE, QFont.Bold))

        subtitle = QLabel("Optimal Speckle Size for DIC Measurements")
        subtitle.setObjectName("subtitle")
        subtitle.setFont(QFont(FONT_FAMILY, FONT_SIZE_SUBTITLE))

        layout.addWidget(title)
        layout.addWidget(subtitle)

        return header

    def create_left_panel(self):
        """Create left panel with inputs"""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setSpacing(15)

        # Camera section
        camera_group = self.create_camera_section()
        layout.addWidget(camera_group)

        # Known parameters section
        known_params_group = self.create_known_params_section()
        layout.addWidget(known_params_group)

        # Calculate button
        self.calc_button = QPushButton("Calculate Speckle Size")
        self.calc_button.setObjectName("calcButton")
        self.calc_button.setFixedHeight(50)
        self.calc_button.setFont(QFont(FONT_FAMILY, FONT_SIZE_HEADING, QFont.Bold))
        self.calc_button.clicked.connect(self.calculate)
        layout.addWidget(self.calc_button)

        layout.addStretch()

        return panel

    def create_camera_section(self):
        """Create camera parameters section"""
        group = QGroupBox("Camera Parameters")
        group.setFont(QFont(FONT_FAMILY, FONT_SIZE_HEADING, QFont.Bold))

        layout = QGridLayout()
        layout.setSpacing(12)
        layout.setContentsMargins(15, 20, 15, 15)

        row = 0

        # Camera selection
        camera_label = QLabel("Camera:")
        camera_label.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(camera_label, row, 0)

        self.camera_combo = QComboBox()
        self.camera_combo.addItems(list(CAMERA_DATABASE.keys()))
        self.camera_combo.setCurrentText("Allied Vision Alvium 1800 U-2050 (IMX541)")
        self.camera_combo.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        self.camera_combo.currentTextChanged.connect(self.on_camera_changed)
        layout.addWidget(self.camera_combo, row, 1)
        row += 1

        # Resolution horizontal
        res_h_label = QLabel("Resolution (H) [px]:")
        res_h_label.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(res_h_label, row, 0)

        self.res_h_input = QLineEdit("4512")
        self.res_h_input.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(self.res_h_input, row, 1)
        row += 1

        # Resolution vertical
        res_v_label = QLabel("Resolution (V) [px]:")
        res_v_label.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(res_v_label, row, 0)

        self.res_v_input = QLineEdit("3211")
        self.res_v_input.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(self.res_v_input, row, 1)
        row += 1

        # Pixel size
        pixel_size_label = QLabel("Pixel Size (µm):")
        pixel_size_label.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(pixel_size_label, row, 0)

        self.pixel_size_input = QLineEdit("2.74")
        self.pixel_size_input.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(self.pixel_size_input, row, 1)
        row += 1

        # Pixels per speckle
        speckle_px_label = QLabel("Pixels per Speckle:")
        speckle_px_label.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(speckle_px_label, row, 0)

        self.speckle_px_range_min = QLineEdit("3")
        self.speckle_px_range_min.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(self.speckle_px_range_min, row, 1)

        self.speckle_px_range_max = QLineEdit("5")
        self.speckle_px_range_max.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(self.speckle_px_range_max, row + 1, 1)

        hint = QLabel("Min - Max pixels")
        hint.setObjectName("hint")
        hint.setFont(QFont(FONT_FAMILY, FONT_SIZE_SMALL))
        layout.addWidget(hint, row + 2, 1)

        layout.setColumnStretch(1, 1)
        group.setLayout(layout)

        # Load initial camera data
        self.on_camera_changed()

        return group

    def on_camera_changed(self):
        """Update fields when camera selection changes"""
        camera_name = self.camera_combo.currentText()
        camera_data = CAMERA_DATABASE[camera_name]

        self.res_h_input.setText(str(camera_data["resolution_h"]))
        self.res_v_input.setText(str(camera_data["resolution_v"]))
        self.pixel_size_input.setText(str(camera_data["pixel_size_um"]))

        # Enable/disable editing for custom camera
        is_custom = camera_name == "Custom"
        self.res_h_input.setReadOnly(not is_custom)
        self.res_v_input.setReadOnly(not is_custom)
        self.pixel_size_input.setReadOnly(not is_custom)

    def create_known_params_section(self):
        """Create known optical parameters section"""
        group = QGroupBox("Known Optical Parameters (Select 2 of 3)")
        group.setFont(QFont(FONT_FAMILY, FONT_SIZE_HEADING, QFont.Bold))

        layout = QGridLayout()
        layout.setSpacing(12)
        layout.setContentsMargins(15, 20, 15, 15)

        row = 0

        # Focal length
        self.focal_check = QCheckBox("Lens Focal Length (mm):")
        self.focal_check.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        self.focal_check.setChecked(True)
        layout.addWidget(self.focal_check, row, 0)

        self.focal_input = QLineEdit("50")
        self.focal_input.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(self.focal_input, row, 1, 1, 2)
        row += 1

        # Object distance
        self.distance_check = QCheckBox("Object Distance (mm):")
        self.distance_check.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        self.distance_check.setChecked(True)
        layout.addWidget(self.distance_check, row, 0)

        self.distance_input = QLineEdit("300")
        self.distance_input.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(self.distance_input, row, 1, 1, 2)
        row += 1

        # Field of view
        self.fov_check = QCheckBox("Field of View (mm):")
        self.fov_check.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        self.fov_check.setChecked(False)
        layout.addWidget(self.fov_check, row, 0)

        self.fov_width_input = QLineEdit("")
        self.fov_width_input.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        self.fov_width_input.setPlaceholderText("Width (mm)")
        layout.addWidget(self.fov_width_input, row, 1)

        layout.setColumnStretch(1, 1)
        layout.setColumnStretch(2, 1)
        group.setLayout(layout)

        # Connect checkboxes
        self.focal_check.stateChanged.connect(self.update_input_states)
        self.distance_check.stateChanged.connect(self.update_input_states)
        self.fov_check.stateChanged.connect(self.update_input_states)

        self.update_input_states()

        return group

    def create_right_panel(self):
        """Create right panel with results"""
        panel = QWidget()
        layout = QVBoxLayout(panel)

        results_group = QGroupBox("Calculation Results")
        results_group.setFont(QFont(FONT_FAMILY, FONT_SIZE_HEADING, QFont.Bold))

        results_layout = QVBoxLayout()
        results_layout.setContentsMargins(15, 20, 15, 15)

        self.results_text = QTextEdit()
        self.results_text.setFont(QFont(FONT_FAMILY_MONO, FONT_SIZE_MONO))
        self.results_text.setReadOnly(True)

        welcome_msg = """
╔════════════════════════════════════════════════════════════╗
║        DIC Speckle Pattern Calculator                      ║
╚════════════════════════════════════════════════════════════╝

Instructions:
  1. Select camera from database or use Custom
  2. Verify resolution (H × V pixels) and pixel size
  3. Set desired pixels per speckle range (typically 3-5)
  4. Select 2 of 3 optical parameters:
     • Lens focal length (mm)
     • Object distance (mm)
     • Field of view width (mm)
  5. Click "Calculate Speckle Size"

The calculator will determine:
  • Missing optical parameter
  • Sensor dimensions (calculated from resolution)
  • Optimal speckle size on specimen
  • Recommended painting techniques

Ready to calculate! 🚀
"""
        self.results_text.setText(welcome_msg)

        results_layout.addWidget(self.results_text)
        results_group.setLayout(results_layout)

        layout.addWidget(results_group)

        return panel

    def update_input_states(self):
        """Update input field states based on checkboxes"""
        checked = [self.focal_check.isChecked(),
                   self.distance_check.isChecked(),
                   self.fov_check.isChecked()]
        num_checked = sum(checked)

        # Enable/disable inputs based on checkbox
        self.focal_input.setEnabled(self.focal_check.isChecked())
        self.distance_input.setEnabled(self.distance_check.isChecked())
        self.fov_width_input.setEnabled(self.fov_check.isChecked())

        # If exactly 2 are checked, disable the other checkboxes
        if num_checked >= 2:
            if not self.focal_check.isChecked():
                self.focal_check.setEnabled(False)
            if not self.distance_check.isChecked():
                self.distance_check.setEnabled(False)
            if not self.fov_check.isChecked():
                self.fov_check.setEnabled(False)
        else:
            self.focal_check.setEnabled(True)
            self.distance_check.setEnabled(True)
            self.fov_check.setEnabled(True)

    def calculate(self):
        """Perform calculation using DICSpeckleCalculator"""
        try:
            # Get camera parameters
            res_width = int(self.res_h_input.text())
            res_height = int(self.res_v_input.text())
            pixel_size_um = float(self.pixel_size_input.text())
            speckle_px_min = int(self.speckle_px_range_min.text())
            speckle_px_max = int(self.speckle_px_range_max.text())

            # Validate selections
            focal_known = self.focal_check.isChecked()
            distance_known = self.distance_check.isChecked()
            fov_known = self.fov_check.isChecked()

            num_known = sum([focal_known, distance_known, fov_known])

            if num_known != 2:
                raise ValueError("Please select exactly 2 of 3 optical parameters")

            # Prepare parameters for calculator
            focal_length_mm = float(self.focal_input.text()) if focal_known else None
            working_distance_mm = float(self.distance_input.text()) if distance_known else None
            fov_width_mm = float(self.fov_width_input.text()) if fov_known else None

            # Use the DICSpeckleCalculator class from the module
            self.calculator = DICSpeckleCalculator(
                resolution_width_px=res_width,
                resolution_height_px=res_height,
                pixel_pitch_um=pixel_size_um,
                pixels_per_speckle_range=(speckle_px_min, speckle_px_max),
                focal_length_mm=focal_length_mm,
                working_distance_mm=working_distance_mm,
                fov_width_mm=fov_width_mm
            )

            # Format results
            results = self._format_results()
            self.results_text.setText(results)

        except ValueError as e:
            error_msg = f"""
╔════════════════════════════════════════════════════════════╗
║                   ⚠️  ERROR                                ║
╚════════════════════════════════════════════════════════════╝

{str(e)}

Please check your inputs:
  • Select exactly 2 of 3 optical parameters
  • All values must be positive numbers
  • Use decimal point (.) for decimals
  • For FOV, provide width value
"""
            self.results_text.setText(error_msg)

        except Exception as e:
            error_msg = f"""
╔════════════════════════════════════════════════════════════╗
║                   ⚠️  ERROR                                ║
╚════════════════════════════════════════════════════════════╝

{str(e)}
"""
            self.results_text.setText(error_msg)

    def _format_results(self):
        """Format calculation results"""
        if not self.calculator:
            return "No calculation performed"

        calc = self.calculator
        wd = calc._calculate_working_distances()

        results = f"""
╔════════════════════════════════════════════════════════════╗
║              CALCULATION RESULTS                           ║
╚════════════════════════════════════════════════════════════╝

📷 CAMERA SETUP
─────────────────────────────────────────────────────────────
  Sensor:             {calc.sensor_w_mm_native:.1f} × {calc.sensor_h_mm_native:.1f} mm
  Resolution:         {calc.resolution_width_px} (H) × {calc.resolution_height_px} (V) px
  Pixel size:         {calc.pixel_pitch_um} µm

📏 OPTICAL PARAMETERS
─────────────────────────────────────────────────────────────
  Focal length:       {calc.focal_length_mm:.3f} mm
  Object distance:    {calc.working_distance_mm:.3f} mm
  Field of view:      {calc.fov_width_mm:.3f} × {calc.fov_height_mm:.3f} mm
  Magnification:      {calc.magnification:.6f}×

  Working distances (optical models):
    Thin Lens:        {wd['thin_lens']:.3f} mm
    Thick Lens (±2%): {wd['thick_lens']:.3f} mm
    BFD (65% of f):   {wd['bfd_65']:.3f} mm
    BFD (70% of f):   {wd['bfd_70']:.3f} mm
    BFD (75% of f):   {wd['bfd_75']:.3f} mm

📊 IMAGE SCALE
─────────────────────────────────────────────────────────────
  Pixel on object:    {calc.pixel_size_mm * 1000:.2f} µm
  Image scale:        {calc.image_scale_px_per_mm:.2f} px/mm

🎯 OPTIMAL SPECKLE SIZE
─────────────────────────────────────────────────────────────
  Range:              {calc.speckle_size_min_mm:.3f} – {calc.speckle_size_max_mm:.3f} mm
  Target size:        {calc.speckle_size_target_mm:.3f} mm ({calc.speckle_size_target_mm * 1000:.1f} µm)
  Recommended subset: {calc.subset_size_min_px} – {calc.subset_size_max_px} px

🎨 RECOMMENDED PAINTING TECHNIQUES
─────────────────────────────────────────────────────────────
"""

        for i, technique in enumerate(calc.recommended_techniques, 1):
            results += f"  {i}. {technique}\n"

        results += """
💡 TIPS
─────────────────────────────────────────────────────────────
  • Use random speckle distribution
  • Aim for 50% black, 50% white coverage
  • Avoid clusters or regular patterns
  • Test correlation quality before full test
"""
        return results

    def apply_stylesheet(self):
        """Apply stylesheet"""
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
        QLineEdit:read-only {{
            background-color: #252525;
            color: {COLOR_TEXT_SECONDARY};
        }}
        QLineEdit:disabled {{
            background-color: #1a1a1a;
            color: {COLOR_TEXT_MUTED};
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
        QCheckBox {{
            spacing: 8px;
            color: {COLOR_TEXT_PRIMARY};
        }}
        QCheckBox::indicator {{
            width: 18px;
            height: 18px;
            border: 2px solid {COLOR_BORDER};
            border-radius: 4px;
            background-color: {COLOR_BG_INPUT};
        }}
        QCheckBox::indicator:checked {{
            background-color: {COLOR_ACCENT};
            border-color: {COLOR_ACCENT};
        }}
        QCheckBox::indicator:disabled {{
            background-color: #1a1a1a;
            border-color: #2a2a2a;
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
        """
        self.setStyleSheet(stylesheet)


def main():
    app = QApplication(sys.argv)
    window = DICSpeckleCalculatorGUI()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()