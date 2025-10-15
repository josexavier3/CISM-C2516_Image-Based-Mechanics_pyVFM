# CISM C2516: José Xavier

Welcome to the repository for the **Virtual Fields Method (VFM)** hands-on session materials from the CISM Advanced School on [Image-Based Mechanics: An Overview of Experimental and Numerical Approaches](https://cism.it/en/activities/courses/C2516/).

## 📚 Course Information

**CISM Advanced School**  
*"Image-Based Mechanics: An Overview of Experimental and Numerical Approaches"*

- **Location:** Udine, Italy
- **Dates:** October 6–10, 2025
- **Course Code:** C2516
- **Coordinators:** Julien Réthoré and José Xavier
- **Lecture:** "The virtual fields method: extracting material parameters from heterogeneous fields: hands-on session"

More information: [https://cism.it/en/activities/courses/C2516/](https://cism.it/en/activities/courses/C2516/)

## 👨‍🔬 Author

**José Xavier**  
Universidade NOVA de Lisboa, NOVA FCT, UNIDEMI  
[https://userweb.fct.unl.pt/~jmc.xavier/](https://userweb.fct.unl.pt/~jmc.xavier/)

## 📖 Theoretical Background

This implementation is based on the book:

*"The Virtual Fields Method: Extracting Constitutive Mechanical Parameters from Full-Field Deformation Measurements"*  
by F. Pierron and M. Grédiac

## 📁 Repository Structure
```
VFM_Unnotched_Iosipescu_Test/
│
├── 1_manuallyVFs/              # Manual virtual fields implementation
│   └── [Implementation files for manually defined VFs]
│
├── 2_piecewiseVFs/             # Piecewise virtual fields approach
│   └── [Implementation files for piecewise VFs]
│
├── ANSYS_APDL/                 # ANSYS APDL finite element models
│   ├── Mtransf3D.mac           # Material transformation macro
│   ├── inputopt.dat            # Input file (23 design variables)
│   └── [Other APDL scripts]
│
└── FEM_Iosipescu_Orthotropic.py  # Main Python implementation
```

### Root Directory Files

- **C2516_Flyer_Rethore_Xavier.pdf** - Course flyer with detailed information
- **Lec06_Xavier.pdf** - Lecture slides
- **Xavier_Introduction.pdf** - Introduction to the VFM methodology
- **.gitignore** - Git ignore configuration

## Installation and Usage

The code in this repository uses Jupyter Notebooks (.ipynb) that can be run through a Jupyter Notebook environment. We recommend using the Jupyter Notebook environment provided by Anaconda.

**To run the code:**
1. Open your Jupyter environment
2. Navigate to the location of the downloaded files
3. Open any of the Python notebooks (.ipynb files) found in the folders
4. The subdirectories contain Python files and experimental data files that will be called by the main notebooks

## Requirements

- Python 3.x
- Jupyter Notebook
- NumPy
- Matplotlib
- (Additional dependencies as specified in the notebooks)

## Reporting Bugs or Issues

Feel free to contact [jmc.xavier@fct.unl.pt] should you encounter any mistakes or bugs in the code!

## Additional Resources

For more information about the CISM course and related materials, visit:  
https://cism.it/en/activities/courses/C2516/


## 📝 License

This code is free for non-profit academic and research use.
It is released under the GNU General Public License v3.0 (GPLv3).
You may use, modify, and share it, as long as you keep the same GPLv3 licence.