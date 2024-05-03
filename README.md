# VP_prediction

VP_prediction is a GUI-based application designed to estimate the saturation vapor pressure of a given molecule. Users input the SMILES notation of a molecule, and the program automatically detects functional groups critical for estimating vapor pressure. It also calculates and displays other relevant molecular parameters, such as the number of carbon (C), hydrogen (H), oxygen (O) atoms, and molecular weight (MW).

# Program Layout

<img width="655" alt="image" src="https://github.com/cyshen93/-VP_prediction/assets/129934466/80182a13-64d9-437a-84f1-3cd2ef2a1d2d">

# Installation and Usage
## Prerequisites
Python v3.0 or higher.
Additional Python modules: tkinter, collections, rdkit, datetime, numpy.

## Setup
Download all the files into a single folder on your computer.
Run the main file named SIMPOL.py.
The main GUI will appear automatically upon execution.

# Features
## Manual Mode
Similar to the original SIMPOL program, you manually input molecular information into the input box, including the number of C, H, O atoms, MW, and the number of different functional groups.Go into this mode by clicking on the manual mode in the main window.

## Auto Mode
When you click the auto-mode, the window will look like the following:
<img width="768" alt="image" src="https://github.com/cyshen93/-VP_prediction/assets/129934466/7e0f4a00-98f2-40d5-970a-d4af15d6e806">

Single Mode: Input the SMILES and temperature in the upper panel, then click "Compute Single Mode". The vapor pressure, C* (saturation concentration), and molecule-related information will be updated in the lower panel. The molecular structure will be displayed on the upper right. If the input SMILES is invalid, a message box will notify you.

Batch Mode: Prepare a text file containing a series of SMILES. Enter the temperature and click "Or Load Your SMILES File" to select your SMILES text file. After clicking "Compute Batch Mode", the program processes all the SMILES in your input file. A message box will appear once processing is complete. The output file is saved in the current folder, named in the format smileYYYYMMDDHHMM.txt, where the numbers represent the save time.


# Acknowledgements
The method to calculate vapor pressure is based on the SIMPOL method proposed by Pankow and Asher (2008). The corresponding paper is included in this folder. This program builds upon the work of Chenyang Shi, which can be found at https://github.com/curieshicy/JRgui.

