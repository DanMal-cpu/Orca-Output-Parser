
# -*- coding: utf-8 -*-
"""
@author: dmalcolm
"""




# ============ GENERAL CHEMISTRY ============

p_block_elements = ["B", "C", "N", "O", "F", "Ne",
                    "Al", "Si", "P", "S", "Cl", "Ar", 
                    "Ga", "Ge", "As", "Se", "Br", "Kr", 
                    "In", "Sn", "Sb", "Te", "I", "Xe", 
                    "Tl", "Pb", "Bi", "Po", "At", "Rn"]

d_block_elements = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
                    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
                    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
                    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn"]



# ============ STRUCTURE DEFINITIONS ============

# adjust and update bond cutoffs as needed, note that sometimes a computed bond may be slightly longer than the benchmarked length value.

bond_cutoffs = {"C": 1.800, "Fe": 3.000, "H": 1.100, "Mo": 3.000, "O": 2.600, "P": 1.725, "W": 3.000}

class_groups = {
    "FRAMEWORK METAL - OVERALL": {"prefix": "all_M", "classifications": ["TM-Bridging Framework Metal", "Regular Framework Metal"]},
    "LIGAND - OVERALL": {"prefix": "lig", "classifications": ["Ligand"]},
    "OXYGEN - OVERALL": {"prefix": "all_O", "classifications": ["Terminal Oxygen", "TM-Terminal Oxygen", "M-M Bridging Oxygen", "M-TM Bridging Oxygen", "M-Het Bridging Oxygen", "TM-Het Bridging Oxygen"]},

    "FRAMEWORK METAL (M-O-M)": {"prefix": "reg_M", "classifications": ["Regular Framework Metal"]},
    "FRAMEWORK METAL (M-O-TM)": {"prefix": "mbridtm_M", "classifications": ["TM-Bridging Framework Metal"]},
    
    "HETEROATOM": {"prefix": "het", "classifications": ["Heteroatom"]},
    
    "LIGAND - OXYGEN (Lig-O-TM)": {"prefix": "tmtolig_O", "classifications": ["TM-Ligand Oxygen"]},  

    "OXYGEN - TERMINAL (M=O)": {"prefix": "term_O", "classifications": ["Terminal Oxygen"]},
    "OXYGEN - TERMINAL (TM=O)": {"prefix": "termTM_O", "classifications": ["TM-Terminal Oxygen"]},
    "OXYGEN - BRIDGING (M-O-M)": {"prefix": "mbridm_O", "classifications": ["M-M Bridging Oxygen"]},
    "OXYGEN - BRIDGING (M-O-TM)": {"prefix": "mbridtm_O", "classifications": ["M-TM Bridging Oxygen"]},
    "OXYGEN - BRIDGING (X-O-M)": {"prefix": "hetbridm_O", "classifications": ["M-Het Bridging Oxygen"]}, 
    "OXYGEN - BRIDGING (X-O-TM)": {"prefix": "hetbridtm_O", "classifications": ["TM-Het Bridging Oxygen"]},
    
    "TM-SUBSTITUTE": {"prefix": "tmsub", "classifications": ["TM-Substitute"]}}



# ============ SYSTEM PATHS ============
    
# orca_dir needs to give the file path to the folder containing orca.exe
# multiwfn_dir needs to do the same but for Multiwfn.exe

orca_dir = r"C:\Users\dmalcolm\ORCA_6.1.0" 
multiwfn_dir = r"C:\Users\dmalcolm\Downloads\Multiwfn_3.8_dev_bin_Win64" 



# ============ TD-DFT PROCESSING OPTIONS ============

# spectral_range can select specific regions of the spectra using "uv", "vis", or "all"
# transition_categories can likewise select groups of transitions based on their nature via "high", "local", or "all"

spectral_range = "all"
transition_categories = "all"
