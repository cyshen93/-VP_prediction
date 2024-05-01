# VP_prediction

This GUI program is to estimate the saturatio vapor pressure of a given molecule. Given the SMILE of a molecule, it can automatically detect the functional groups used to estimate the vpor pressure. The other parameters like the number of C, H, O and moleweight will also be calculated and displayed.

The program layout will look like the following:

<img width="655" alt="image" src="https://github.com/cyshen93/-VP_prediction/assets/129934466/80182a13-64d9-437a-84f1-3cd2ef2a1d2d">

# How to use it?

You need to install python v3.0+ in your computer and the modules imported in this program include: tkinter, collections, rdkit, datetime and numpy. Download all the files into one folder in your computer and run the main file called "SIMPOL.py". The main GUI will automatically appear.

In the main GUI, you have two options to calculate the vapor pressure of a molecule or molecules: manual-mode and auto-mode. The manual-mode is similar to the SIMPOL igor program and you need to input all the molecular information into the box like C, H, O, MW, number of different functional groups. The auto-mode is a new feature compared to the old SIMPOL igor program. The only information you need to provide is the SMILE of the molecule. The program can detect all required information itself!

When you click the auto-mode, the window will look like the following:
<img width="768" alt="image" src="https://github.com/cyshen93/-VP_prediction/assets/129934466/7e0f4a00-98f2-40d5-970a-d4af15d6e806">
The auto-mode also has two options: single-mode and batch-mode. Single-mode is straightforward and you need to input the SMILES, Temp in the upper panel, and then click on the "compute single-mode". The vapor pressure, Cstar and molecule-related information will be updated on the lower panel. The molecule structure will also be displayed on the upper right. If the input SMILE is not good to use, a message box will be prompted to let you know.

For the Batch-mode, you need to prepare a text file which contains a series of SMILES. Then you need to input the Temp and click on the "Or Load Your SMILES File" button to select the SMILES .txt file. You need to click on the "compute batch-mode" button and the program will process all the smiles in your input file. When it's done, it will prompt a message box to notify you. Keep in mind that the output file is saved in the current folder with a name like "smile202404292100.txt'. The numbers in the name correspond to the save time. 

# Acknowledgements

The SIMPOL method to calculate the vapor pressure is proposed by Pankow and Asher, 2008 and the corresponding paper is attached in this folder. Part of my code is build on Chenyang Shi's work : https://github.com/curieshicy/JRgui.


