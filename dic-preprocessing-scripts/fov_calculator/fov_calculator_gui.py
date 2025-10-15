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
    QGridLayout, QSplitter, QFrame
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
FONT_MONO_FALLBACK = "JetBrains Mono, Consolas, Courier New, monospace"

FONT_SIZE_TITLE = 22
FONT_SIZE_SUBTITLE = 11
FONT_SIZE_HEADING = 12
FONT_SIZE_NORMAL = 10
FONT_SIZE_SMALL = 9
FONT_SIZE_MONO = 10

# ============================================================================
# GLOBAL COLOR CONFIGURATION - Grey Button Variant
# ============================================================================
COLOR_ACCENT = "#F59E0B"           # ← Change this for orange header/borders
COLOR_ACCENT_HOVER = "#FBBF24"
COLOR_ACCENT_PRESSED = "#D97706"
COLOR_ACCENT_LIGHT = "#FCD34D"

# Button colors - GREY BASED ⭐
COLOR_BUTTON_BG = "#3F3F46"        # Change from #52525B to #3F3F46
COLOR_VIZ_FOV = "#9CA3AF"          # Change from #F59E0B to #9CA3AF
COLOR_TEXT_PRIMARY = "#ffffff"     # Change from #E5E7EB to #ffffff

COLOR_BG_MAIN = "#1a1a1a"          # ← Main background
COLOR_BG_SECONDARY = "#242424"     # ← Secondary panels
COLOR_BG_INPUT = "#2d2d2d"         # ← Input field background
COLOR_BG_HEADER = "#1f1f1f"        # ← Header background

COLOR_BORDER = "#404040"           # ← Borders
COLOR_BORDER_FOCUS = "#F59E0B"     # ← Focused border
COLOR_SEPARATOR = "#404040"        # ← Separator line

COLOR_TEXT_PRIMARY = "#E5E7EB"     # ← Main text color
COLOR_TEXT_SECONDARY = "#9CA3AF"   # ← Secondary text
COLOR_TEXT_MUTED = "#6B7280"       # ← Hint text

COLOR_VIZ_FOV = "#F59E0B"          # ← FOV rectangle color (orange)
COLOR_VIZ_ROI = "#06B6D4"          # ← ROI rectangle color (cyan)

# ============================================================================


class FOVVisualizer(FigureCanvas):
    """Matplotlib canvas for FOV visualization"""

    def __init__(self, parent=None):
        self.fig = Figure(figsize=(6, 4), facecolor=COLOR_BG_SECONDARY)
        super().__init__(self.fig)
        self.setParent(parent)
        self.ax = self.fig.add_subplot(111, facecolor=COLOR_BG_SECONDARY)
        self.setup_plot()

    def setup_plot(self):
        """Initialize the plot"""
        self.ax.set_aspect('equal')
        self.ax.grid(True, alpha=0.3, color='white', linestyle='--')
        self.ax.set_xlabel('Width (mm)', color='white', fontsize=FONT_SIZE_NORMAL)
        self.ax.set_ylabel('Height (mm)', color='white', fontsize=FONT_SIZE_NORMAL)
        self.ax.tick_params(colors='white', labelsize=FONT_SIZE_SMALL)

        for spine in self.ax.spines.values():
            spine.set_color('white')

        self.fig.tight_layout()

    def plot_fov(self, fov_width, fov_height, roi_width, roi_height):
        """Plot FOV and ROI rectangles"""
        self.ax.clear()
        self.setup_plot()

        # Plot FOV rectangle (orange)
        fov_rect = plt.Rectangle(
            (-fov_width / 2, -fov_height / 2), fov_width, fov_height,
            fill=False, edgecolor=COLOR_VIZ_FOV, linewidth=2.5,
            label='Field of View'
        )
        self.ax.add_patch(fov_rect)

        # Plot ROI rectangle (green)
        roi_rect = plt.Rectangle(
            (-roi_width / 2, -roi_height / 2), roi_width, roi_height,
            fill=True, facecolor=COLOR_VIZ_ROI, alpha=0.15,
            edgecolor=COLOR_VIZ_ROI, linewidth=2, linestyle='--',
            label='Recommended ROI (95%)'
        )
        self.ax.add_patch(roi_rect)

        # Add center crosshair
        self.ax.plot([0, 0], [-fov_height / 2 * 1.1, fov_height / 2 * 1.1],
                     'w--', alpha=0.3, linewidth=0.8)
        self.ax.plot([-fov_width / 2 * 1.1, fov_width / 2 * 1.1], [0, 0],
                     'w--', alpha=0.3, linewidth=0.8)

        # Set limits with padding
        padding = 1.15
        self.ax.set_xlim(-fov_width / 2 * padding, fov_width / 2 * padding)
        self.ax.set_ylim(-fov_height / 2 * padding, fov_height / 2 * padding)

        # Add title
        self.ax.set_title('Field of View Visualization', color='white',
                          fontsize=FONT_SIZE_HEADING, fontweight='bold', pad=10)

        # Add legend OUTSIDE the plot area (top right)
        legend = self.ax.legend(
            loc='upper left',
            bbox_to_anchor=(1.02, 1),
            facecolor=COLOR_BG_SECONDARY,
            edgecolor='white',
            fontsize=FONT_SIZE_SMALL,
            labelcolor='white',
            framealpha=0.9
        )
        for text in legend.get_texts():
            text.set_fontfamily(FONT_FAMILY)

        # Adjust layout to prevent legend cutoff
        self.fig.tight_layout()

        self.draw()

class FOVCalculatorGUI(QMainWindow):

    def __init__(self):
        super().__init__()

        self.sensor_dimensions = {
            "1/3.2\" (4.5×3.4mm)": (4.54, 3.42),
            "1/2.5\" (5.8×4.3mm)": (5.76, 4.29),
            "1/2.3\" (6.2×4.6mm)": (6.17, 4.55),
            "2/3\" (8.8×6.6mm)": (8.8, 6.6),
            "1.0\" (12.8×9.6mm)": (12.8, 9.6),
            "1.1\" (13.2×8.8mm)": (13.2, 8.8),
            "4/3\" (17.3×13.0mm)": (17.3, 13.0),
        }

        self.load_custom_fonts()
        self.init_ui()
        self.apply_stylesheet()

    def load_custom_fonts(self):
        """Load custom fonts from system or local directory"""
        font_dir = Path(__file__).parent / "fonts"
        if font_dir.exists():
            for font_file in font_dir.glob("*.ttf"):
                QFontDatabase.addApplicationFont(str(font_file))
            for font_file in font_dir.glob("*.otf"):
                QFontDatabase.addApplicationFont(str(font_file))

        available_fonts = QFontDatabase.families()
        if FONT_FAMILY not in available_fonts:
            print(f"⚠️  Font '{FONT_FAMILY}' not found. Using fallback.")
        if FONT_FAMILY_MONO not in available_fonts:
            print(f"⚠️  Font '{FONT_FAMILY_MONO}' not found. Using fallback.")

    def init_ui(self):
        """Initialize the user interface"""
        self.setWindowTitle("FOV Calculator - DIC Setup")
        self.setGeometry(100, 100, 1200, 700)

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

        # Right panel - Results Text
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

        title = QLabel("🔬 Field of View Calculator")
        title.setObjectName("title")
        title.setFont(QFont(FONT_FAMILY, FONT_SIZE_TITLE, QFont.Bold))

        subtitle = QLabel("Digital Image Correlation Setup Tool")
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
        self.calc_button = QPushButton("Calculate FOV")
        self.calc_button.setObjectName("calcButton")
        self.calc_button.setFixedHeight(50)
        self.calc_button.setFont(QFont(FONT_FAMILY, FONT_SIZE_HEADING, QFont.Bold))
        self.calc_button.clicked.connect(self.calculate_fov)
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

        # Sensor size
        sensor_label = QLabel("Sensor Size:")
        sensor_label.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(sensor_label, 0, 0)

        self.sensor_combo = QComboBox()
        self.sensor_combo.addItems(list(self.sensor_dimensions.keys()))
        self.sensor_combo.setCurrentText("1.1\" (13.2×8.8mm)")
        self.sensor_combo.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(self.sensor_combo, 0, 1)

        # Focal length
        focal_label = QLabel("Focal Length (mm):")
        focal_label.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(focal_label, 1, 0)

        self.focal_input = QLineEdit("50")
        self.focal_input.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        self.focal_input.setPlaceholderText("e.g., 50")
        layout.addWidget(self.focal_input, 1, 1)

        focal_hint = QLabel("Typical: 16, 25, 35, 50, 75, 100")
        focal_hint.setObjectName("hint")
        focal_hint.setFont(QFont(FONT_FAMILY, FONT_SIZE_SMALL))
        layout.addWidget(focal_hint, 2, 1)

        # Stand-off distance
        standoff_label = QLabel("Stand-off Distance (mm):")
        standoff_label.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        layout.addWidget(standoff_label, 3, 0)

        self.standoff_input = QLineEdit("500")
        self.standoff_input.setFont(QFont(FONT_FAMILY, FONT_SIZE_NORMAL))
        self.standoff_input.setPlaceholderText("e.g., 500")
        layout.addWidget(self.standoff_input, 3, 1)

        standoff_hint = QLabel("Distance from lens to specimen")
        standoff_hint.setObjectName("hint")
        standoff_hint.setFont(QFont(FONT_FAMILY, FONT_SIZE_SMALL))
        layout.addWidget(standoff_hint, 4, 1)

        layout.setColumnStretch(1, 1)
        group.setLayout(layout)

        return group

    def create_visualization_section(self):
        """Create visualization display"""
        group = QGroupBox("Field of View Visualization")
        group.setFont(QFont(FONT_FAMILY, FONT_SIZE_HEADING, QFont.Bold))

        layout = QVBoxLayout()
        layout.setContentsMargins(10, 20, 10, 10)

        self.fov_canvas = FOVVisualizer()
        layout.addWidget(self.fov_canvas)

        group.setLayout(layout)
        return group

    def create_right_panel(self):
        """Create right panel with text results"""
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
║        Welcome to FOV Calculator                       ║
╚════════════════════════════════════════════════════════╝

Instructions:
  1. Select sensor size
  2. Enter focal length (mm)
  3. Enter stand-off distance (mm)
  4. Click "Calculate FOV"

Ready to calculate! 🚀
"""
        self.results_text.setText(welcome_msg)

        layout.addWidget(self.results_text)
        group.setLayout(layout)

        return group

    def calculate_fov(self):
        """Perform FOV calculation"""
        try:
            sensor_str = self.sensor_combo.currentText()
            focal_length = float(self.focal_input.text())
            standoff_distance = float(self.standoff_input.text())

            if focal_length <= 0 or standoff_distance <= 0:
                raise ValueError("Values must be positive")

            if standoff_distance <= focal_length:
                raise ValueError("Stand-off distance must be greater than focal length")

            L_CS_width, L_CS_height = self.sensor_dimensions[sensor_str]

            # Calculate FOV
            L_FOV_width = L_CS_width * ((standoff_distance - focal_length) / focal_length)
            L_FOV_height = L_CS_height * ((standoff_distance - focal_length) / focal_length)
            L_FOV_diagonal = np.sqrt(L_FOV_width ** 2 + L_FOV_height ** 2)
            magnification = L_CS_width / L_FOV_width

            # Format results
            results = f"""
╔════════════════════════════════════════════════════════╗
║              CALCULATION RESULTS                       ║
╚════════════════════════════════════════════════════════╝

📥 INPUT PARAMETERS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Sensor:             {sensor_str}
  Sensor dimensions:  {L_CS_width:.2f} × {L_CS_height:.2f} mm
  Focal length:       {focal_length:.1f} mm
  Stand-off distance: {standoff_distance:.1f} mm

📊 FIELD OF VIEW
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  FOV Width:          {L_FOV_width:.2f} mm
  FOV Height:         {L_FOV_height:.2f} mm
  FOV Diagonal:       {L_FOV_diagonal:.2f} mm
  Magnification:      {magnification:.4f}×

🎯 RECOMMENDED ROI (90-95% of FOV)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Width range:        {L_FOV_width * 0.90:.2f} - {L_FOV_width * 0.95:.2f} mm
  Height range:       {L_FOV_height * 0.90:.2f} - {L_FOV_height * 0.95:.2f} mm

💡 TIPS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  • Higher magnification = smaller FOV
  • Specimen should fit in 90-95% of FOV
  • Consider lighting setup requirements
"""

            self.results_text.setText(results)

            # Update visualization
            roi_width = L_FOV_width * 0.95
            roi_height = L_FOV_height * 0.95
            self.fov_canvas.plot_fov(L_FOV_width, L_FOV_height, roi_width, roi_height)

        except ValueError as e:
            error_msg = f"""
╔════════════════════════════════════════════════════════╗
║                   ⚠️  ERROR                            ║
╚════════════════════════════════════════════════════════╝

{str(e)}

Please check your inputs:
  • Values must be positive numbers
  • Stand-off must be > focal length
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
            background-color: {COLOR_ACCENT};
            border: none;
            border-radius: 8px;
            color: {COLOR_TEXT_PRIMARY};
            padding: 10px;
            font-family: {FONT_FAMILY}, {FONT_FALLBACK};
        }}
        QPushButton#calcButton:hover {{
            background-color: {COLOR_ACCENT_HOVER};
        }}
        QPushButton#calcButton:pressed {{
            background-color: {COLOR_ACCENT_PRESSED};
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
    window = FOVCalculatorGUI()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()