'''
        ttk.Label(self.bottom_frame, text = 'Carbon Groups', font = ('Arial', 16, 'bold')).grid(row = 0, column = 2, padx = 5, pady = 5)
        
        ttk.Label(self.bottom_frame, text = 'non-aromatic rings',font = ('Arial',16, 'bold')).grid(row = 1, column =2, padx = 5, pady = 5)
        ttk.Entry(self.bottom_frame, width = 5, state = 'readonly',textvariable = self.G[1][0]).grid(row = 1, column =3, padx = 5, pady = 5)
       
        ttk.Label(self.bottom_frame, text = 'aromatic rings', font = ('Arial', 16, 'bold')).grid(row = 2, column = 2, padx = 5, pady = 5)
        ttk.Entry(self
        '''

"""
SIMPOL in python version
This code is written by Chuanyang Shen (cyshen2012@gmail.com) Apr 29th 2024
This script cross-reference SIMPOL in igor version (Pankow and Asher, 2008 ACP.), jrgui.py by Chenyang shi and Joback_Hand_Pick.py 
"""

import os
import tkinter.ttk as ttk
import tkinter as tk
from tkinter import *
from rdkit import Chem
from rdkit.Chem import Descriptors
from tkinter import filedialog
from collections import Counter
from operator import itemgetter
import numpy as np
from datetime import datetime
from rdkit.Chem import Draw

my_color = 'grey91'
my_size = 14
# my_color = "#34A2FE"

class Overall_Look:
    def __init__(self, master):
        self.master = master
        self.master.title('SIMPOL GUI')
        self.master.configure(background = my_color)
        self.master.minsize(300, 400) # width + height
        self.master.resizable(False, False)

        self.style = ttk.Style()
        self.style.configure('TFrame', background = my_color)
        self.style.configure('TButton', background = my_color)
        self.style.configure('TLabel', background = my_color, font = ('Arial', 16))
        self.style.configure('Header.TLabel', font = ('Arial', 24, 'bold'))
        ## the top frame
        self.frame_header = ttk.Frame(master, relief = RIDGE, padding = (30, 15))
        self.frame_header.pack()
        ttk.Label(self.frame_header, text = 'Welcome to SIMPOL GUI program!', style = 'Header.TLabel').grid(row = 0, column =2, columnspan = 4)
        ttk.Label(self.frame_header, wraplength = 600, text = ('This is the PYTHON GUI for SIMPOL. This program is used to predict the saturation vapor pressure and Cstar of a given molecular structure. This method is proposed by Pankow and Asher, 2008. This program runs in two modes -- manual and auto modes. In manual mode, the users need to input the number of different functional groups of a molecule. In auto mode, all parameters can be detected automatically if a SMILES string of the molecule is given. Auto mode also includes batch mode in which a txt file containing a series of SMILES can be given and the predictions will be saved into a txt ouput file.')).grid(row = 1, column =2, columnspan =4)

        ## the bottom frame
        self.bottom_frame = ttk.Frame(master, relief = FLAT, padding = (30,15))
        self.bottom_frame.pack()

        # two buttons to choose from
        Button(self.bottom_frame, text = 'Manual Mode', relief = RAISED, command = self.hand_pick, font = ('Arial', 20, 'bold'),foreground = 'blue').grid(row=0, column =0, columnspan=3, padx = 5, pady = 5)
        Button(self.bottom_frame, text = 'Auto Mode', relief = RAISED, command = self.smiles, font = ('Arial', 20, 'bold'),fg = 'blue').grid(row=0, column = 3, columnspan =3, padx = 5, pady=5)

    def hand_pick(self):
        self.hand_pick = tk.Toplevel(self.master)
        self.GUI = Manual_Mode(self.hand_pick)

    def smiles(self):
        self.smiles = tk.Toplevel(self.master)
        self.GUI = Auto_Mode(self.smiles)

class Manual_Mode:
    def __init__(self, master):
        self.master = master


class Auto_Mode:
    def __init__(self, master):
        self.master = master
        self.master.title('SIMPOL Auto-Mode')
        self.master.configure(background = my_color)
        self.master.minsize(700, 700) # width + height
        self.master.resizable(False, True)


        self.style = ttk.Style()
        self.style.configure('TFrame', background = my_color)
        self.style.configure('TButton', background = my_color)
        self.style.configure('TCheckbutton', background = my_color)
        self.style.configure('TRadiobutton', background = my_color)
        self.style.configure('Tlabel', background = my_color, font = ('Arial', 16))
        self.style.configure('Header.TLabel', font = ('Arial', 24, 'bold'))

        ## load image
        __location__ = os.path.realpath(os.path.join(os.getcwd(),os.path.dirname(__file__)))
        self.logo = PhotoImage(file = os.path.join(__location__, 'structure0.gif')).subsample(1,1)


        self.frame_header = ttk.Frame(master, relief = RIDGE, padding = (30, 15))
        self.frame_header.pack()

        ## top frame to put Entry box for SMILES and Temperature
        self.top_frame = ttk.Frame(self.master, padding = (30, 15),relief = FLAT)
        self.top_frame.pack()

        self.lbl_img = ttk.Label(self.top_frame, image = self.logo)
        self.lbl_img.grid(row = 0, column =2, rowspan = 7,sticky ='e')
       
        lbl_smile = ttk.Label(self.top_frame, text = 'Input your SMILES', font = ('Arial', 16, 'bold')).grid(row =0, column =0, padx = 5, pady = 5, )
        self.SMILES = StringVar()
        self.SMILES.set('CC1=CC(=C(C=C1)C(=O)O)O')
        ent_smile = ttk.Entry(self.top_frame, width = 30, font = ('Arial', 14), textvariable = self.SMILES).grid(row = 1, column = 0,padx = 5, pady = 5)
        btn_smile = ttk.Button(self.top_frame, text = 'Or Load Your SMILES File', command = self.popoutwindow, style = 'TButton').grid(row = 2, column = 0, padx = 5, pady = 5)


        lbl_temp = ttk.Label(self.top_frame, text = 'Temperature (K)',justify = RIGHT, font = ('Arial', 16, 'bold')).grid(row = 3, column = 0,padx = 5, pady = 5)
        self.Temp = StringVar()
        self.Temp.set(298)
        ent_temp = ttk.Entry(self.top_frame, width = 20, font = ('Arial', 14), textvariable = self.Temp).grid(row = 4, column = 0, padx = 5, pady = 5)

        ## middle frame to put the Calculation Button
        # self.middle_frame = ttk.Frame(self.master, padding = (30, 15))
        # self.middle_frame.pack()
        Button(self.top_frame, text = 'Compute single-mode', command = self.calculate_vp_single, relief = RAISED, font = ('Arial',15,'bold'),foreground = 'blue').grid(row =5, column = 0, padx = 5) 
        Button(self.top_frame, text = 'Compute batch-mode', command = self.calculate_vp_batch,relief = RAISED, font = ('Arial',15,'bold'),foreground = 'blue').grid(row = 6, column = 0,  padx = 5)

        ## bottom frame to put the calcualtion results
        self.bottom_frame = ttk.Frame(self.master,width = 200, height = 110, relief = SUNKEN)
        self.bottom_frame.pack(padx = 20, pady = 10)
       
               
        # tk.Label(self.bottom_frame, text = 'Results', font = ('Arial', my_size, 'bold'), fg = 'blue').grid(row = 0, column = 0, padx = 5, pady = 5)
        ttk.Label(self.bottom_frame, text = 'VPatm', font = ('Arial',my_size, 'bold')).grid(row = 1, column =0, padx = 5, pady = 5)
        self.lbl_vp = tk.Label(self.bottom_frame, text = 'None', font = ('Arial', my_size, 'bold'),bg = 'white',width =8 )
        self.lbl_vp.grid(row = 1, column = 1, columnspan = 1, padx = 5, pady = 5)
        
        ttk.Label(self.bottom_frame, text = 'Cstar', font = ('Arial',my_size, 'bold')).grid(row = 2, column =0, padx = 5, pady = 5)
        # ttk.Entry(self.bottom_frame, width = 5, state = 'readonly',textvariable = self.Cstar).grid(row = 2, column = 1, padx = 5, pady = 5)
        self.lbl_Cstar = tk.Label(self.bottom_frame, text = 'None', font = ('Arial', my_size, 'bold'), bg = 'white',width = 8)
        self.lbl_Cstar.grid(row =2, column = 1, columnspan = 1, padx = 5, pady =5) 

        ttk.Label(self.bottom_frame, text = 'VP (Evaporation)', font = ('Arial', my_size, 'bold')).grid(row = 3, column = 0 , padx = 5, pady = 5)
        self.lbl_vp_evap = tk.Label(self.bottom_frame, text = 'None', font = ('Arial', my_size, 'bold'), bg = 'white', width = 8)
        self.lbl_vp_evap.grid(row = 3, column = 1, columnspan = 1, padx = 5, pady = 5)

        ttk.Label(self.bottom_frame, text = 'Cstar (Evaporation)', font = ('Arial', my_size, 'bold')).grid(row = 4, column =0, padx = 5, pady = 5)
        self.lbl_cstar_evap = tk.Label(self.bottom_frame, text = 'None', font = ('Arial', my_size, 'bold'), bg = 'white', width = 8)
        self.lbl_cstar_evap.grid(row = 4, column = 1, padx = 5, pady = 5)
        ## put the formula and ring results
        
       
        # tk.Label(self.bottom_frame, text = 'Formula', font = ('Arial', my_size, 'bold'), fg = 'blue').grid(row = 4, column = 0, padx = 5, pady = 5)
        label = ['C', 'H', 'O', 'MW']
                
        for i in range(4):
            ttk.Label(self.bottom_frame, text = label[i], font = ('Arial',my_size, 'bold')).grid(row = i+6, column = 0, padx = 5, pady =  5)
        
        self.lbl_numC = tk.Label(self.bottom_frame, text = '0', font = ('Arial', my_size, 'bold'),bg = 'white', width = 6 )
        self.lbl_numC.grid(row = 6, column = 1, columnspan = 1, padx = 5, pady = 5)
        self.lbl_numH = tk.Label(self.bottom_frame, text = '0', font = ('Arial', my_size, 'bold'),bg = 'white',width = 6 )
        self.lbl_numH.grid(row = 7, column = 1, columnspan = 1, padx = 5, pady = 5)
        self.lbl_numO = tk.Label(self.bottom_frame, text = '0', font = ('Arial',my_size, 'bold'),bg = 'white',width =6 )
        self.lbl_numO.grid(row = 8, column = 1, columnspan = 1, padx = 5, pady = 5)
        self.lbl_MW = tk.Label(self.bottom_frame, text = '0', font = ('Arial', my_size, 'bold'), bg = 'white', width = 6)
        self.lbl_MW.grid(row = 9, column = 1, columnspan = 1, padx = 5, pady = 5)
        
        # tk.Label(self.bottom_frame, text = 'Carbon Groups', font = ('Arial', my_size, 'bold'), fg = 'blue').grid(row = 10, column = 0, padx = 5, pady = 5) 

        label = ['non-aromatic rings', 'aromatic rings', 'non-aromatic C=C', 'aromatic C=C']
        name = ['carb1','carb2', 'carb3','carb4']
       
        for i in range(4):
            ttk.Label(self.bottom_frame, text = label[i], font = ('Arial', my_size, 'bold')).grid(row = i+10, column = 0, padx = 5, pady = 5)

        self.lbl_carb1 = tk.Label(self.bottom_frame, text = '0', font = ('Arial', my_size, 'bold'),bg = 'white', width = 6 )
        self.lbl_carb1.grid(row = 10, column = 1, columnspan = 1, padx = 5, pady = 5)
        self.lbl_carb2 = tk.Label(self.bottom_frame, text = '0', font = ('Arial', my_size, 'bold'),bg = 'white',width = 6 )
        self.lbl_carb2.grid(row = 11, column = 1, columnspan = 1, padx = 5, pady = 5)
        self.lbl_carb3 = tk.Label(self.bottom_frame, text = '0', font = ('Arial', my_size, 'bold'),bg = 'white',width = 6 )
        self.lbl_carb3.grid(row = 12, column = 1, columnspan = 1, padx = 5, pady = 5)
        self.lbl_carb4 = tk.Label(self.bottom_frame, text = '0', font = ('Arial', my_size, 'bold'), bg = 'white', width = 6)
        self.lbl_carb4.grid(row = 13, column = 1, columnspan = 1, padx = 5, pady = 5)

        ## put the oxygen groups
        # tk.Label(self.bottom_frame, text = 'Oxygen Groups', font = ('Arial', my_size, 'bold'), fg = 'blue').grid(row = 0, column = 3, padx = 5, pady = 5)
        label = ['-OH (alkyl)', '-OH (aromatic)', '-C(O)H (aldehyde)', '-C(O)- (ketone)', '-C(O)OH (acid)', '-COOH (hydroperoxide)','-C(O)OOH (carbonylperoxy acid)', '-C(O)O- (esther)', '-O- (alkyl ether)', '-O- (alicyclic ether)', '-O- (aromatic ether)', '-OO- (peroxide)', '-ONO2 (nitrate)', '-NO2 (nitro)']
        for i in range(14): 
            ttk.Label(self.bottom_frame, text = label[i], font = ('Arial', my_size, 'bold')).grid(row = i, column = 3, padx = 5, pady = 5)
      
        self.lbl_oxy = [None]*14
        for i in range(14):
            self.lbl_oxy[i] = tk.Label(self.bottom_frame, text = '0', font = ('Arial', my_size, 'bold'), bg = 'white', width = 6)
            self.lbl_oxy[i].grid(row = i, column = 4, padx = 8, pady = 6)

                
    def search_functional_group(self, smile):
        """
        Counts the number of different functional groups for SIMPOL calculation in a molecule.

        Parameters:
        - smiles (str): The SMILES string of the molecule.

        Returns:
        - tuple: A tuple containing the number of Formula, Carbon Groups, and Oxygen Group, respectively.
        - The tuple has the form of (G1, G2, G3)
    
        list G1: size of 4, representing the number of C, H, O, moleweigth
        list G2: size of 4, representing the number of non-aromatic rings, aromatic rings, non-aromatic C=C, aromatic C=C
        list G3: size of 14, representing the number of functional groups 
        (0)-OH (alkyl),
        (1) -OH (aromatic),
        (2) -C(O)H (aldehyde), 
        (3) -C(O)- (ketone), 
        (4) -C(O)OH (acid), 
        (5) -COOH (hydroperoxid), 
        (6) -C(O)OOH (carbonylperoxy acid) 
        (7) -C(O)O- (ester), 
        (8) -O- (alkyl ether), 
        (9) -O- (alicyclic ether), 
        (10) -O- (aromatic ether), 
        (11) -OO- (peroxide), 
        (12) -ONO2 (nitrate), 
        (13) -NO2 (nitro)

        """
        # print(smile)

        # Convert the SMILES string to an RDKit molecule object
        m = Chem.MolFromSmiles(smile)
        m2 = Chem.AddHs(m)
        vals = Descriptors.CalcMolDescriptors(m)

        # Part I: calculate B1
        # C, H, O, molweight
        G1 = [0] * 4 # initialize the list
        atom_counts = {atom: 0 for atom in ['C','H','O']}
        # Iterate over all atoms in the molecule
        for atom in m2.GetAtoms():
            # Get the atom symbol
            symbol = atom.GetSymbol()
            # If the symbol is in our list, increment its count in the dictionary
            if symbol in atom_counts:
                atom_counts[symbol] += 1

        G1[0] = atom_counts['C'] 
        G1[1] = atom_counts['H'] 
        G1[2] = atom_counts['O'] 
        G1[3] = Descriptors.MolWt(m2)



        # part II: calculate G2
        # non-aromatic rings; aromatic rings; C=C (non-aromatic); C=C (aromatic)
        G2 = [0]*4  # initialize the list 

        non_aromatic_double_bonds = 0
        aromatic_double_bonds = 0

        # Iterate over all bonds in the molecule
        for bond in m.GetBonds():
            # Check if the bond is a double bond
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetBeginAtom().GetSymbol() == 'C' and bond.GetEndAtom().GetSymbol() == 'C':
                # Check if the bond is aromatic
                if bond.GetIsAromatic():
                    aromatic_double_bonds += 1
                else:
                    non_aromatic_double_bonds += 1
    
        G2[1] = vals['NumAromaticRings']
        G2[0] = vals['RingCount'] - G2[1]
        G2[2] = non_aromatic_double_bonds 
        G2[3] = aromatic_double_bonds 
    

        # part III: calculate B3
        # oxygen groups: 14
    
        # find the smart structure for different functional groups
        smarts = ["[$([CX2H0](=*)=*)]", "[$([CX2H1]#[!#7])]", "[$([CX2H0]#[!#7])]", "[OX2H]-[C]=O", "[#6X3H0;!$([#6X3H0](~O)(~O)(~O))](=[#8X1])[#8X2H0]", "[$([#6X3H0](=[OX1]));!$([#6X3](=[#8X1])~[#8X2]);R]=O", "[CH;D2;$(C-!@C)](=O)", "[OX2H;!$([OX2H]-[#6]=[O]);!$([OX2H]-a)]", "[O;H1;$(O-!@c)]", "[#8X2H0;R;!$([#8X2H0]~[#6]=[#8])]","[$([CX3H0](=[OX1]));!$([CX3](=[OX1])-[OX2]);!R]=O", "[OX2H0;!R;!$([OX2H0]-[#6]=[#8])]", "[$([#7X3,#7X3+][!#8])](=[O])~[O-]", "[OX1H0;!$([OX1H0]~[#6X3]);!$([OX1H0]~[#7X3]~[#8])]", "[#7X2H0;R]", "[#7X3H1;R]", "[#7X2H1]","[#7X2H0;!R]","[#6X2]#[#7X1H0]","[NX3H2]", "[NX3H1;!R]", "[#7X3H0;!$([#7](~O)~O)]","[SX2H]","[#16X2H0;!R]","[#16X2H0;R]", "[R;CX3H1,cX3H1]", "[$([R;#6X3H0]);!$([R;#6X3H0]=[#8])]","[R;CX4H2]","[R;CX4H]","[R;CX4H0]", "[CX3H2]", "[!R;CX3H1;!$([CX3H1](=O))]","[$([!R;#6X3H0]);!$([!R;#6X3H0]=[#8])]","[CX4H3]","[!R;CX4H2]", "[!R;CX4H]","[!R;CX4H0]","[F]","[Cl]","[Br]","[I]"]

        smarts1 = [0]*14
        smarts1[0] = smarts[7] #'[#6][OX2H]' # smarts[7]  # -OH alkyl
        smarts1[1] = smarts[8] #'[OX2H][cX3]:[c]' #smarts[8]  # -OH phenol
        smarts1[2] = '[CX3H1](=O)[#6]' # smarts[6] denotes alkyl aldehyde  # -C(O)H aldehyde
        smarts1[3] = '[#6][CX3](=O)[#6]'   # ketone
        smarts1[4] = smarts[3]   # -C(O)OH (acid)
        smarts1[5] = '[OX2H][OX2]'#ester/hydroperoxide '[OX2H][OX2]' 
        smarts1[6] ='[CX3](=[OX1])[OX2][OX2H1]' # carbonylperoxy acid

        smarts1[7] = smarts[4]  # -C(O)O- ester
        smarts1[8] = "[OX2H0;!R;!$([OX2H0]-[#6]=[#8])]" #smarts[11] # "[OX2H0;!R;!$(O[N])]" # smarts[11]  # -O- (non-ring/akyl ether)
        smarts1[9] = smarts[9]  #-O- (ring)
        smarts1[10] = '[#6][OX2H0][c]'       # aromatic ether -O-
        smarts1[11] = '[#6][OX2,OX1-][OX2,OX1-][#6]'  # -OO- (peroxide)
        smarts1[12] = '[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]' # -ONO2- (nitrate)
        smarts1[13] = smarts[12] # '([NO3](=[OX1])(=[OX1])(=[OX1])   '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]' # -NO2 (nitro)

        # initialize the parameter
        tuples = []
        index_list = []
        final_index_and_length =[]

        # Iterate over all 14 functional groups
        for index, smart in enumerate(smarts1):
            if m.HasSubstructMatch(Chem.MolFromSmarts(smart)) == True:
                tuples.append(m.GetSubstructMatches(Chem.MolFromSmarts(smart)))
                index_list.append(index)

        temp = self.return_non_duplicate_index(tuples)

        for i in temp:
            final_index_and_length.append([index_list[i[0]], i[1]])

        G3 = [0]*14
    

        for i in final_index_and_length:
            index = i[0]
            num = i[1]
            G3[index] = num

        G3[0] = vals['fr_Al_OH'] # alkyl -OH
        G3[1] = vals['fr_Ar_OH'] # aromatic -OH
        G3[2] = vals['fr_aldehyde']
        G3[3] = vals['fr_ketone'] # ketone
        G3[13] = vals['fr_nitro'] # -NO2 nitro group

        return (G1,G2,G3)

        

    def calculate_vp_single(self):
        smile = self.SMILES.get()
        smile = str(smile)
        Temp = self.Temp.get()
        Temp = int(Temp)

        try:
            m = Chem.MolFromSmiles(smile)
            img = Draw.MolToImage(m)
        except:
            tk.messagebox.showwarning("Warning!", "It seems like your SMILES cannot be properly parsed, please check!")

        # m = Chem.MolFromSmiles(smile)
        # img = Draw.MolToImage(m)
        img.save('structure.gif')
        __location__ = os.path.realpath(os.path.join(os.getcwd(),os.path.dirname(__file__)))
        self.logo = PhotoImage(file = os.path.join(__location__, 'structure.gif')).subsample(1,1)

        # self.logo = PhotoImage(file = 'structure.gif').subsample(1,1)
        self.lbl_img['image'] = self.logo
        
        VPatm, Cstar, G = self.calculate_vp(smile, Temp)

        self.VPatm = VPatm
        self.Cstar = Cstar
        self.G = G

        self.lbl_vp['text'] = f"{VPatm:.2e}"
        self.lbl_Cstar['text'] = f"{Cstar:.2e}"

        self.lbl_numC['text'] = f"{self.G[0][0]}"
        self.lbl_numH['text'] = f"{self.G[0][1]}"
        self.lbl_numO['text'] = f"{self.G[0][2]}"
        self.lbl_MW['text'] = f"{self.G[0][3]:.0f}"
        self.lbl_carb1['text'] = f"{self.G[1][0]}"
        self.lbl_carb2['text'] = f"{self.G[1][1]}"
        self.lbl_carb3['text'] = f"{self.G[1][2]}"
        self.lbl_carb4['text'] = f"{self.G[1][3]}"
        for i in range(14):
            self.lbl_oxy[i]['text'] = f"{self.G[2][i]}"



    def calculate_vp(self, smile, Temp):

        """
        Compute the saturation vapor pressure of a molecule.

        Parameters:
        - smiles (str): The SMILES string of the molecule.
        - Temp (K): The temperature (K).

        Returns:  
        - A tuple containing vp and Cstar: (vp, Cstar)
        - vp: unit atm, saturation vapor pressure 
        - Cstar: unit ug/m3 
        """

        # calculate the number of functional groups
        G1,G2,G3 =  self.search_functional_group(smile)

        # B values from Pankow and Asher, 2008 ACP.
        B1 = [None]*31
        B1[0:11] = [-426.938,-411.248,-146.442,35.0262,-87.277,5.73335,-261.268,-725.373,-729.501,-13.7456,-798.796]
        B1[11:21] = [-393.345,-144.334,40.5265,-70.7406,-783.648,-563.872,-453.961,37.1375,-503.710,-35.9763]
        B1[21:31] = [-609.432,-102.367,-1938.02,-5.26919,-284.042,150.093,-20.3387,-838.064,-52.7934,-1615.20]

        B2 = [None]*31
        B2[0:11] = [0.289223,0.896919,1.54528,-0.920839,1.78059,0.0169764,-0.763282,0.826326,0.986017,0.523486,-1.09436]
        B2[11:21] = [-0.951778,-1.85617,-2.4378,-1.06674,-1.03439,-0.718416,-0.326105,-2.66753,1.04092,-0.408458]
        B2[21:31] = [1.50436,-0.716253,0.648262,0.306435,-0.625424,0.0239875,-5.48718,-1.096,-0.463689,0.901669]

        B3 = [None]*31
        B3[0:11] = [0.00442057,-0.00248607,0.00171021,0.00224399,-0.00307187,-0.000628957,-0.00168213,0.00250957,-0.00292664,0.000550298,0.00524132]
        B3[11:21] = [-0.00219071,-2.37481E-05,0.00360133,0.00373104,-0.00107148,0.00263016,-0.00013978,1.01483E-3,-4.12746E-3,1.67264E-3]
        B3[21:31] = [-9.09024E-4,-2.9067E-4,1.73245E-3,3.25397E-3,-8.22474E-4,-0.00337969,0.00839075,-0.000424385,-5.11647E-3,1.44536E-3]

        B4 = [None] * 31
        B4[0:11] = [0.292846,0.140312,-0.278291,-0.09363,-0.104341,0.00755434,0.289038,-0.232304,0.178077,-0.27695,-0.22804]
        B4[11:21] = [0.305843,0.28829,0.0986422,-0.144003,0.315535,-0.049947,-0.0393916,0.214233,0.18279,-0.0998919]
        B4[21:31] = [-0.135495,-0.588556,0.034794,-0.681506,-0.080549,0.0152769,0.107884,0.281812,0.384965,0.266889]

        B1 = np.array(B1)
        B2 = np.array(B2)
        B3 = np.array(B3)
        B4 = np.array(B4)
    
        Tref = 293.15   
        Bref = B1/Tref + B2 + B3*Tref + B4*np.log(Tref)
        BT = B1/Temp + B2 + B3*Temp + B4*np.log(Temp)

        Gcontrib = [0] * 31
        Gcontrib[0] = 1
        #                 C,    g2,    g3       g4      g5,      g6,   g7,    g8,    g9,   g10 
        Gcontrib[1:11] = [G1[0],0,     G2[1],  G2[0],  G2[2],   G2[3], G3[0], G3[2], G3[3],G3[4]]               
        Gcontrib[11:21] =[G3[7],G3[8], G3[9],  G3[10], G3[12],  G3[13],G3[1], 0,      0,    0   ]
        Gcontrib[21:31] =[0,    0,     0,      0,      0,       G3[11],G3[5], G3[6],  0, 0   ] 
        Gcontrib = np.array(Gcontrib)

        VP = sum(BT*Gcontrib)
        VPatm = 10**(VP)
        VPtorr= VPatm*760

        Cstar = (VPtorr *1000000*G1[3])/(760*Temp*8.21e-5)
        
        G = [G1,G2,G3]

        return (VPatm, Cstar,G)



    def calculate_vp_batch(self):
        smiles = self.smile_string
        Temp = self.Temp.get()
        Temp = float(Temp)

        num = len(smiles)
        vp = [None]*num
        cstar = [None]*num
        G = [None]*num
        for i, smile in enumerate(smiles):
            if not Chem.MolFromSmiles(smile):
                continue
            else:
                vp[i], cstar[i], G[i] = self.calculate_vp(smile, Temp)

        self.VPatm = vp
        self.Cstar = cstar 
        self.G = G
        
        ct = datetime.now()
        fn = f"smile_{ct.year:4d}{ct.month:02d}{ct.day:02d}_{ct.hour:02d}{ct.minute:02d}.txt"
        # fn = 'smile.txt'
        f = open(fn, "a")
        f.write("smile  Temp   Vp(atm)      Cstar   Formula     Carbon_Groups   oxygen_groups\n")
        for i in range(num):
            if vp[i]:
                tline = f"{smiles[i]}   {Temp}  {vp[i]} {cstar[i]}  {G[i][0][:]}    {G[i][1][:]}    {G[i][2][:]} \n" 
            else:
                tline = f"{smiles[i]}   None    None    None    None    None    None \n"
            f.write(tline)

        f.close()

        tk.messagebox.showinfo(title = None, message = 'Output txt file has been successfully saved to the current path!')

        return None


    def popoutwindow(self):
        self.input_file_name = filedialog.askopenfile(defaultextension = '.txt', mode = 'r', filetypes = [('Text Documents', '*.txt')])
        self.text = self.input_file_name.readlines()
        print(self.text)
        self.smile_string = [] # a list to store the SMILES strings

        for line in self.text:
            self.smile_string.append(line.strip('\n'))
       
        print(self.smile_string)

    def return_non_duplicate_index(self, tuples):
        # given a list of sets return index of non_duplicate items

        ## step 1, create a new tuple, named "new_tuples"
        new_tuples = []
        for i in tuples:
            for j in i:
                new_tuples.append(set(j))


        ## step 2, creat a dictionary storing one to one relationship between new_tuple and old_tuple
        values = []
        for index, item in enumerate(tuples):
            if len(item) == 1:
                values.append(index)
            else:
                for i in [index]*len(item):
                    values.append(i)

        keys = [i for i in range(len(new_tuples))]
        dict_tuples = {}
        for i, j in zip(keys, values):
            dict_tuples[i] = j

        ## step 3: remove duplicates in sets terminology
        remove_index = []
        for index_1, item in enumerate(new_tuples):
            for index_2 in range(index_1 + 1, len(new_tuples)):
                if len(item & new_tuples[index_2]) != 0:
                    if len(item) > len(new_tuples[index_2]):
                        remove_index.append(index_2)
                    elif len(item) < len(new_tuples[index_2]):
                        remove_index.append(index_1)
                    elif len(item) == len(new_tuples[index_2]):
                        remove_index.append(index_2)
        remain_sets = set(range(len(new_tuples))).difference(set(remove_index))


        ## step 4: spit out final index and length
        index_1 = []
        index_length = []

        for i in remain_sets:
            index_1.append(dict_tuples[i])

        ## count
        counts = Counter(index_1)
        list_counts = counts.most_common()
        for i in range(len(list_counts)):
            index_length.append([list_counts[i][0], list_counts[i][1]])

        index_length = sorted(index_length, key = itemgetter(0))
        return index_length



def main():
    # create a new window
    root = Tk()

    # creat an instance of Auto_Mode class
    GUI = Overall_Look(root)
    # GUI = Auto_Mode(root)

    # event loop: listens for events like click/input/... and blocks any code that comes after it from running until you close the window
    root.mainloop()

if __name__ == '__main__': main()






