
# -*- coding: utf-8 -*-
"""
@author: dmalcolm
"""

import os
       
import Config as C     
import Electronic_Properties as EP
import ORCA_Commands as OC
import Physical_Properties as PP      
import TDDFT_Multiwfn_Comb as TMC  

from glob import glob
from shutil import move

from Molecule_Classification import classify_atoms, count_elements, desig_prop
from Parameter_Storage import OrcaParameterStorage
from TDDFT_Transmittance import generate_spectrum

parameter_storage = OrcaParameterStorage()
    
#Fe3+_Mo_Keggin_H-t-butyl-Peroxo.out
#Mo_H-t.txt




def calc_analysis(filename):
    """ Extract all data from ORCA output file. Returns dict of results. """    
    
    results = {}
    
    parameter_storage.store_parameter("filename", filename)                                                             # Entire filename is made available to OrcaParameterStorage() class.
    
    
    # ============ BASIC STRUCTURE DATA ============
    
    results["coords"] = PP.extract_and_label_coordinates(filename)                                                      # Retrieves optimized set of coordinates from file.
    parameter_storage.store_parameter("coords", results["coords"])
    
    atom_data = PP.store_atom_coordinates(results["coords"])                                                            # Splits coordinates into separate x, y, and z dictionaries for easier retrieval, analysis and sorting.
    results["atom_id_list"], results["x_coords"], results["y_coords"], results["z_coords"] = atom_data
    
    
    # ============ BONDING DATA ============
    
    results["bond_lengths"] = PP.calculate_bond_lengths(results["coords"], C.bond_cutoffs)                              # Stores interatomic distances based on maximum lengths allowed by bond_cutoffs, "creates" bonds.
    parameter_storage.store_parameter("bond_lengths", results["bond_lengths"])  
    
    results["bond_orders"] = PP.extract_bond_orders(filename)                                                           # Retrieves Mayer bond orders directly from output file.
    parameter_storage.store_parameter("bond_orders", results["bond_orders"])
    
    
    # ============ ELECTRONIC PROPERTIES ============
   
    results["charge"] = EP.extract_mol_charge(filename)                                                                 # Finds molecular charge from input section of output file.
    parameter_storage.store_parameter("charge", results["charge"])
    
    results["qa_charges"] = EP.extract_qa_charges(filename, list(results["x_coords"].keys()))                           # Retrieves Mulliken Gross Atomic Charges (QA) directly from output file.
    parameter_storage.store_parameter("qa_charges", results["qa_charges"])
    
    results["sort_qa"] = EP.sort_dict_ntop(results["qa_charges"])                                                       # Sorts QA Charges based on charge (-ve to +ve).
    parameter_storage.store_parameter("sort_qa", results["sort_qa"])
    parameter_storage.store_parameter("sorted_qa_all", results["sort_qa"])
    
    results["hir_charges"] = EP.extract_hir_charges(filename, list(results["x_coords"].keys()))                         # Retrieves Hirshfeld Charges directly from output file.
    parameter_storage.store_parameter("hir_charges", results["hir_charges"])
    
    results["sort_hir"] = EP.sort_dict_ntop(results["hir_charges"])                                                     # Sorts Hirshfeld Charges (-ve to +ve).
    parameter_storage.store_parameter("sort_hir", results["sort_hir"])
    
    results["hir_spin_pop"] = EP.extract_hir_pop(filename, list(results["x_coords"].keys()))                            # Retrieves Hirschfeld Spin values directly from output file.
    parameter_storage.store_parameter("hir_spin_pop", results["hir_spin_pop"])
    
    results["sort_hir_spin_pop"] = EP.sort_dict_pton(results["hir_spin_pop"])                                           # Sorts Hirschfeld Spin values (+ve to -ve).
    parameter_storage.store_parameter("sort_hir_spin_pop", results["sort_hir_spin_pop"])
    
    results["loew_spin_pop"] = EP.extract_loew_pop(filename, list(results["x_coords"].keys()))                          # Retrieves Löwdin Spin values directly from output file.
    parameter_storage.store_parameter("loew_spin_pop", results["loew_spin_pop"])
    
    results["sort_loew_spin_pop"] = EP.sort_dict_pton(results["loew_spin_pop"])                                         # Sorts Löwdin Spin values (+ve to -ve).
    parameter_storage.store_parameter("sort_loew_spin_pop", results["sort_loew_spin_pop"])
    
    results["loew_orbital_spin"] = EP.extract_orb_spin_contribution(filename, list(results["x_coords"].keys()))         # Retrieves Reduced Orbital Spin Populations directly from output file. Collects values from valence orbitals only (p for C, d for Fe...). 
    parameter_storage.store_parameter("loew_orbital_spin", results["loew_orbital_spin"])
    
    orbital_spin_values = {atom: data["total_orbital_spin"] for atom, data in results["loew_orbital_spin"].items()} 
    results["sort_loew_orbital_spin"] = EP.sort_dict_pton(orbital_spin_values)                                          # Sorts Total Reduced Orbital Spin Population values (+ve to -ve).
    parameter_storage.store_parameter("sort_loew_orbital_spin", results["sort_loew_orbital_spin"])
    
    results["frontier"] = EP.extract_homo_lumo_join(filename)                                                           # Retrieves HOMO and LUMO orbital energies and designations. Accounts for spin polarisation.
    parameter_storage.store_parameter("frontier", results["frontier"])
    
    results["extended_frontier"] = EP.extended_extract_homo_lumo_join(filename)                                         # Retrieves HOMO-10 to LUMO+10 orbital energies and designations. Accounts for spin polarisation.
    parameter_storage.store_parameter("extended_frontier", results["extended_frontier"])
    
    
    # ============ CLASSIFICATION ============
    
    results["element_list"] = count_elements(results["atom_id_list"])                                                   # Counts numbers of each atom by element.
    results["atom_class"] = classify_atoms(results["element_list"], results["bond_lengths"])                            # Classify all atoms based on element and number/ types of bonds.
    parameter_storage.store_parameter("atom_class", results["atom_class"])
    
    results = desig_prop(results, parameter_storage)                                                                    # Collects all data for a specific atom type.
    
    
    # ============ ENERGY DATA ============
    
    results["energy_data"] = EP.extract_energy_values(filename)                                                         # Retrieves Electronic Energy and Thermal Correction values, stores Calculation Type.
    parameter_storage.store_parameter("energy_data", results["energy_data"])
    
    if results["energy_data"]["calculation_type"] in ["Frequency", "Numerical Frequency"]:
        results["vib_freq"] = EP.extract_vib_freq(filename, results["energy_data"]["calculation_type"])                 # If Frequency calculation, retrieves Vibrational Frequencies. Reports any Imaginary Frequencies.
        parameter_storage.store_parameter("vib_freq", results["vib_freq"])

    return results



def multiwfn_nto_parsing():
    """ Process TD-DFT data and generate NTO and CDD cube files for relevant transitions using Multiwfn. """
    
    code_dir = os.getcwd()                                                                                          # Directory containing code.
 
    matches = []
    molden_files = []
    
    
    # ============ FIND ORCA FILES IN SPYDER DIRECTORY ============
    
    print(f"\n{'='*75}")
    print("📊 Step 1: Scanning for ORCA files in Spyder directory ...")
    print(f"{'='*75}\n")
                     
    with TMC.working_directory(code_dir):
        for gbw_file in [f for f in os.listdir() if f.endswith(".gbw")]:                                            # Find all .gbw files.
            base_name = gbw_file[:-4]                                                                               # Save filename without the .gbw suffix.

            if os.path.exists(f"{base_name}.out"):                                                                  # Check to see if a .out file with that name exists in the same folder.
                matches.append({"base_name": base_name})                                                            # Store filename without extension to "matches" for later use.
                print(f"   ✅ Found matching pair: {base_name}.out / {base_name}.gbw")
            else:
                print(f"   ⚠️  No matching .out for: {gbw_file}")
    
    if not matches:
        print("No valid ORCA files found. Exiting.")
        return
    
        
    # ============ GENERATE MOLDEN FILES ============
    
    print(f"\n\n{'='*75}")
    print("📊 Step 2: Copying files to ORCA directory and generating Molden files...")
    print(f"{'='*75}\n")
    

    for match in matches:
        base_name = match["base_name"]                                                                              # Access filenames from "matches".
        
        if OC.copy_orca_files_to_orca_dir(code_dir, C.orca_dir, base_name):                                         # Copy .gbw and .out file pairs to from code to orca directory.
            if OC.generate_molden_file(C.orca_dir, base_name):                                                      # Generate .molden file per pair.
                molden_files.append(base_name)                                                                      # Store filename of each new .molden file.
            else:
                print(f"   ❌ Failed to generate molden file for {base_name}")              
        else:
            print(f"   ❌ Failed to copy files for {base_name}")
    
    if not molden_files:
        print("No files were successfully copied and processed!")
        return
        
    current_dir = os.getcwd()                                                                                       # Create Orbitals folder in the code directory.
    orbitals_base = os.path.join(current_dir, "Orbitals")
    os.makedirs(orbitals_base, exist_ok=True)
    
    for filename in molden_files:                                                                                   # Process each molden file.
        
        calc_folder = os.path.join(orbitals_base, filename)                                                         # Create the calculation folder in Orbitals.
        os.makedirs(calc_folder, exist_ok=True)
        
        
        # ============ GENERATE TD-DFT UV-VIS SPECTRA ============
        
        print(f"\n\n{'='*75}")
        print("📊 Step 3: Generate UV-Vis Spectra for each TD-DFT File...")
        print(f"{'='*75}\n")
        
        generate_spectrum(f"{filename}.out")                                                                        # Generate UV-Vis spectrum.
            
        pattern = f"{filename}_absorbance_sigma*.png"                                                               # Find and move the matching .png spectra file.
        matching_files = glob(pattern)
    
        if matching_files:
            for png_file in matching_files:
                move(png_file, os.path.join(calc_folder, png_file))
                print(f"   📊 UV-Vis file moved to: {os.path.join(calc_folder, png_file)}")
        else:
            print(f"   ⚠️  No UV-Vis .png file found for {filename}")
    
            
        # ============ PARSE TD-DFT FILE FOR SIGNIFICANT TRANSITIONS ============
        
        print(f"\n\n{'='*75}")
        print("📊 Step 4: Parsing Output File for Transitions...")
        print(f"{'='*75}")
        
        try:
            abs_data = TMC.parse_orca_tddft_absorption(f"{filename}.out")                                           # Collect simulated absorption data from TD-DFT output file.
            print(f"\n   🔍 Found {len(abs_data)} transitions")
            
            if C.spectral_range == "uv":
                abs_data = [t for t in abs_data if 200 <= t["wavelength"] <= 400]
                print(f"\n   🔍 Filtered to UV region: {len(abs_data)} transitions in this range")
            elif C.spectral_range == "vis":
                abs_data = [t for t in abs_data if 400 <= t["wavelength"] <= 800]
                print(f"\n   🔍 Filtered to Vis region: {len(abs_data)} transitions in this range")
            
            selected = TMC.select_excitation_peaks_by_criteria(abs_data, 0.1)                                       # Select peaks from absorption spectra based on specific intensity conditions.
            nto_data = TMC.parse_nto_for_selected_states(f"{filename}.out", selected)                               # Parse NTO data for each selected excitation and collect enough NTO excitations for fair representation.
            TMC.save_selected_excitations(selected, nto_data, filename, calc_folder)                                # Save all collected peaks and NTO data in "Orbitals" folder in .txt file.
            
            transition_numbers = TMC.get_transition_numbers_by_category(selected, C.transition_categories)          # Collect NTO transitions into a dictionary. This function is key for processing selected NTO data with Multiwfn.
            #transition_numbers = [3]                                                                               # !!! Alternatively use this for selecting specific NTOs/ peaks. Useful for quick code tests !!!
            
                
            # ============ GENERATE NTO FILES FOR SIGNIFICANT TRANSITIONS ============
            
            print(f"\n   ➡️Need NTO workup for {len(transition_numbers)} unique transitions: {transition_numbers}")
            
            TMC.run_multiwfn_generation(C.multiwfn_dir, C.orca_dir, filename, transition_numbers, nto_data)         # Generate NTO and CDD cube files using Multiwfn.
            OC.cleanup_temp_files(C.orca_dir, filename)                                                             # Cleanup any temporary, leftover files.
            
            print(f"\n✅ Finished Processing: {filename}\n\n")
            
        except Exception as e:
            print(f"❌ Error parsing TDDFT output: {e}")
            return
    
    
    for file in molden_files:
        OC.cleanup_orca_files(C.orca_dir, file)                                                                     # Remove .out and .gbw files from ORCA directory.
        
    print(f"\n{'='*75}\n🎉 All Files Processed Successfully!\n{'='*75}")



def main():
    print("ORCA Data Parser UI")
    print("-------------------\n")
    
    filename = input("Enter ORCA output filename: ").strip()
    
        
    # ============ PROCESS AS STANDARD OR TD-DFT OUTPUT ============
    
    tddft_present = TMC.check_tddft_present(filename)
    if tddft_present:
        try:
            multiwfn_nto_parsing()
            return
        except Exception as e:
            return
    
    calc_analysis(filename)
    
        
    # ============ ASK FOR USER INPUT ============
    
    while True:
        selection = input("\nEnter:\n - parameter name (see ReadMe or enter 'list')\n - 'save' to export data \n - 'quit'\n ").strip()
        
        if selection.lower() == "quit":
            break
        elif selection.lower() == "list":
            print("\nAvailable parameters:")
            for param in sorted(parameter_storage.parameters.keys()):
                print(f"  {param}")
        elif selection.lower() == "save":
            filename = input("Enter output filename (default: orca_parameters.txt): ").strip()
            if not filename:
                filename = "orca_parameters.txt"
            parameter_storage.save_to_raw_file(filename)                                        # Generates summary txt file containing all raw data, sorted and grouped accordingly.
            trim_name = os.path.splitext(filename)[0]
            parameter_storage.save_to_trimmed_file(trim_name)                                   # Generates more concise txt file containing more key and relevant results.

        else:
            value = parameter_storage.get_parameter(selection)                                  # Allows for viewing specific parameters directly in the console window.
            if value is not None:
                print(f"\n{selection}:")

                formatted_value = parameter_storage.format_parameter(selection, value)
                print(formatted_value)
            else:
                print(f"\nParameter '{selection}' not found.")

if __name__ == "__main__":
    main()