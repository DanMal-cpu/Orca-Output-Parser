# -*- coding: utf-8 -*-
"""
@author: dmalcolm
"""

import os

from contextlib import contextmanager
from re import search
from shutil import move
from subprocess import run, TimeoutExpired
from time import sleep, time




photon_constant = 1239.841984           # eV * nm


def check_tddft_present(filename):
    """ Quickly scan the output file for TDDFT keyword. """
    
    try:
        with open(filename, "r") as f:
            content = f.read()
            if "tddft" in content.lower():
                return True
            return False
    except FileNotFoundError:
        print(f"⚠️ File {filename} not found!")
        return False



def generate_nto_cube_files(multiwfn_dir, orca_dir, trans_num, interpreted_ntos):
    """ Generate hole (H) and electron (L) cube files for a transition using its .mwfn file. Multiwfn generates cube files with generies names that are renamed after generation. """

    multiwfn_exe = None
    with working_directory(multiwfn_dir):                                                                           # Locate Multiwfn executable.
        for exe in ["Multiwfn.exe", "Multiwfn_noGUI.exe", "Multiwfn.bat"]:
            if os.path.exists(exe):
                multiwfn_exe = os.path.join(multiwfn_dir, exe)
                break
    
    if not multiwfn_exe:
        print("   ❌ Multiwfn executable not found for NTO cube generation")
        return False
    
    if not interpreted_ntos:
        print("   ⚠️  No NTO data provided")
        return False
    
    success_count = 0
    total_expected = len(interpreted_ntos) * 2                                                                      # Hole + Electron for each NTO.
    mwfn_full_path = os.path.join(orca_dir, f"S{trans_num}.mwfn")
    
    print(f"\n   🔄 Generating cube files for {len(interpreted_ntos)} NTOs...")

    for nto_idx, nto in enumerate(interpreted_ntos, 1):                                                             # Generate NTO-associated cube file for each electron and hole of a specific NTO.
        print(f"\n   NTO{nto_idx}:")
    
        for desig, suffix in [(nto["electron_desig"], "electron"), (nto["hole_desig"], "hole")]:
            
            # Commands to pass to Multiwfn
            commands = "\n".join([
                "",                         # First ENTER
                "",                         # Second ENTER
                "5",                        # 5 Output and plot specific property within a spatial region (calc. grid data)
                "4",                        # 4 Value of orbital wavefunction 
                desig,                      # Input orbital index e.g. h for HOMO
                "2",                        # 2 Medium quality grid, covering whole system, about 512000 points in total
                "2",                        # 2 Export data to a Gaussian-type cube file in current folder
                "0",                        # 0 Return to main menu
                "q"                         # Quit
            ]) + "\n"
    
            try:
                run([multiwfn_exe, mwfn_full_path], input=commands, capture_output=True, text=True, cwd=orca_dir, timeout=120)
                    
                for file in os.listdir(orca_dir):                                                                   # Find and rename the generated cube file.
                    if file.endswith(".cub") and "CDD" not in file and f"NTO{nto_idx}_{suffix}" not in file:
                        if time() - os.path.getmtime(os.path.join(orca_dir, file)) < 10:                            # Check if file was just created (within last 10 seconds).
                            new_name = f"S{trans_num}_NTO{nto_idx}_{suffix}.cub"
                            os.rename(os.path.join(orca_dir, file), os.path.join(orca_dir, new_name))
                            print(f"      ✅ {new_name} generated")
                            success_count += 1
                            break

            except Exception as e:
                print(f"      ❌ Error generating {suffix} cube: {str(e)}")
    
    if success_count == total_expected:
        print(f"\n   ✅ Generated all {total_expected} cube files for S{trans_num}")
        return True
    else:
        print(f"   ⚠️  Only generated {success_count}/{total_expected} cube files for S{trans_num}")
        return False



def generate_transition_files(multiwfn_dir, orca_dir, base_filename, trans_num, workspace_info, nto_data, idx, len_trans):
    """ Run Multiwfn to generate NTOs and CDD for specified transitions. """
    
    print(f"\n\n{'='*75}")
    print(f"📊 Step 6: Running Multiwfn NTO Generation for S{trans_num}   [Transition {idx}/{len_trans}]")
    print(f"{'='*75}\n")
    
    interpreted_ntos = []
    if trans_num in nto_data and nto_data[trans_num]["nto_transitions"]:
        print(f"   🗒️ NTO Data for S{trans_num}:\n")
        print(f"      Total occupation: {nto_data[trans_num]['total_occupation']:.4f}")
        print(f"      Number of NTOs needed: {nto_data[trans_num]['num_nto_needed']}\n")
        print("      NTO transitions (interpreted):")
        
        data = nto_data[trans_num]
        for occ in data["nto_transitions"]:
            elec, hole = interpret_nto_orbital(occ["from"], occ["to"])
            print(f"      {occ['from']}a -> {occ['to']}a (occ={occ['occupation']:.6f})")
            interpreted_ntos.append({"from": occ["from"], "to": occ["to"], "occupation": occ["occupation"], "electron_desig": elec, "hole_desig": hole})
    else:
        print("   ⚠️  No NTO data available")
        return False
    
    multiwfn_exe = None                                                                                         # Locate Multiwfn executable and molden file.
    with working_directory(multiwfn_dir):
        for exe in ["Multiwfn.exe", "Multiwfn_noGUI.exe", "Multiwfn.bat"]:
            if os.path.exists(exe):
                multiwfn_exe = os.path.join(multiwfn_dir, exe)
                break
    
    if not multiwfn_exe:
        print(f"❌ Error: Multiwfn executable not found in {multiwfn_dir}")
        return False
    
    molden_file = os.path.join(orca_dir, f"{base_filename}.molden")                                             # Check if molden file exists.
    if not os.path.exists(molden_file):
        print(f"❌ Error: Molden file not found: {molden_file}")
        return False
    
    
    commands = "\n".join([                                                                                      # Command list for .mwfn and .cdd generation for this single transition.
        "",                     # First ENTER - dismiss settings.ini warning
        "",                     # Second ENTER - confirm file load
        "18",                   # Electron excitation analysis
        "6",                    # Generate NTOs
        "",                     # ENTER to see list
        str(trans_num),         # Select specific excited state
        "3",                    # Output NTO orbitals to .mwfn file
        f"S{trans_num}.mwfn",   # Filename
        "0",                    # Return to main menu
        "18",                   # Electron excitation analysis
        "1",                    # Analyze hole & electron distribution, transition density, and transition electric/magnetic dipole moment density
        "",                     # ENTER to continue
        str(trans_num),         # Select specific excited state
        "1",                    # Visualize and analyze hole, electron, transition density and so on
        "2",                    # Medium quality grid, covering whole system, about 512000 points in total
        "15",                   # Output cube file of charge density difference to current folder
        "0",                    # Return
        "0",                    # Return
        "0",                    # Return
        "q"                     # Quit
    ]) + "\n"

    
    print(f"\n   🔄 Generating S{trans_num}.mwfn and S{trans_num}-CDD.cub ...\n")
    
    try:
        result = run([multiwfn_exe, molden_file], input=commands, capture_output=True, text=True, cwd=orca_dir, timeout=240)
       
        with open(os.path.join(orca_dir, f"{base_filename}_logfile.log"), "w", encoding="utf-8") as f:          # Save complete output to log file.
            f.write(f"\n=== TRANSITION {trans_num} ===\n{commands}\n{result.stdout}\n{result.stderr}")
            
        sleep(1)
        
        if os.path.exists(os.path.join(orca_dir, f"S{trans_num}.mwfn")):                                        # Check if .mwfn file was created.
            size = os.path.getsize(os.path.join(orca_dir, f"S{trans_num}.mwfn"))
            print(f"   ✅ S{trans_num}.mwfn generated ({size:,} bytes)")
        else:
            print(f"   ❌ S{trans_num}.mwfn not found")
            return False
        
        cdd_path = os.path.join(orca_dir, "CDD.cub")                                                            # Rename CDD.cub immediately.
        if os.path.exists(cdd_path):
            move(cdd_path, os.path.join(orca_dir, f"S{trans_num}-CDD.cub"))
            size = os.path.getsize(os.path.join(orca_dir, f"S{trans_num}-CDD.cub"))
            print(f"   ✅ S{trans_num}-CDD.cub generated ({size:,} bytes)")
        else:
            print(f"\n   ⚠️  No CDD.cub file found for transition {trans_num}")
        
        generate_nto_cube_files(multiwfn_dir, orca_dir, trans_num, interpreted_ntos)                            # Generate H and L cube files.  
        
        print("\n   📁 Organising files...")                                                                   # Organize all files into workspace.
        organize_generated_files(orca_dir, trans_num, workspace_info)
        
        print(f"\n   ✅ Successfully processed S{trans_num}\n")
        return True
    
    except TimeoutExpired:
        print(f"   ❌ Timeout for transition {trans_num}")
        return False
    except Exception as e:
        print(f"   ❌ Exception: {str(e)}")
        import traceback
        traceback.print_exc()
        return False



def get_transition_numbers_by_category(selected, categories="all"):
    """ Get transition numbers filtered by category.
    Return sorted list of those transitions. """
    
    category_map = {"high": "high_intensity", "local": "local_peaks"}
    
    if categories == "all":
        transitions = set()
        for key in category_map.values():
            if key in selected and selected[key]:
                for trans in selected[key]:
                    transitions.add(trans["transition"])
        return sorted(list(transitions))
        
    transitions = set()
    
    if isinstance(categories, list):
        for cat in categories:
            if cat in category_map:
                key = category_map[cat]
                if key in selected and selected[key]:
                    for trans in selected[key]:
                        transitions.add(trans["transition"])
    
    return sorted(list(transitions))



def interpret_nto_orbital(from_orb, to_orb):
    """ Interpret NTO orbital numbers and return the Multiwfn selection type. 
    Returns electron and hole designations (e.g. "h", "l")
    
    NOTE: orbital numbers in Multiwfn are 1 integer off what they should be.
    Referring to them using frontier orbital notation provides more consistent results and cuts down on confusion.
    Just keep in mind NTO orbitals are different than frontier orbitals. """
    
    gap = to_orb - from_orb
    
    if gap == 1:
        return("h", "l")
    elif gap == 3:
        return("h-1", "l+1")
    elif gap == 5:
        return("h-2", "l+2")
    elif gap == 7:
        return("h-3", "l+3")
    
    else:
        return(f"orbital {from_orb}", f"orbital {to_orb}")



def organize_generated_files(orca_dir, trans_num, workspace_info):
    """ Move generated files for a specific transition to organised folders. """
    
    trans_folder = workspace_info["transition_folders"].get(trans_num)
    if not trans_folder:
        print(f"   ❌ No workspace folder for S{trans_num}")
        return False
    
    os.makedirs(trans_folder, exist_ok=True)                                # Create folder if it doesn't exist.
    
    files_moved = []
    
    nto_source = os.path.join(orca_dir, f"S{trans_num}.mwfn")               # Move NTO file.
    nto_dest = os.path.join(trans_folder, f"S{trans_num}.mwfn")
    if os.path.exists(nto_source):
        move(nto_source, nto_dest)
        files_moved.append(f"S{trans_num}.mwfn")
        print(f"   📦 Moved: S{trans_num}.mwfn -> S{trans_num}/")
    
    cdd_source = os.path.join(orca_dir, f"S{trans_num}-CDD.cub")            # Move CDD file.
    cdd_dest = os.path.join(trans_folder, f"S{trans_num}-CDD.cub")
    if os.path.exists(cdd_source):
        move(cdd_source, cdd_dest)
        files_moved.append(f"S{trans_num}-CDD.cub")
        print(f"   📦 Moved: S{trans_num}-CDD.cub -> S{trans_num}/")
        
    nto_cube_patterns = ["_hole.cub", "_electron.cub"]                      # Move all NTO cube files.
    for file in os.listdir(orca_dir):
        if file.startswith(f"S{trans_num}_NTO") and any(pattern in file for pattern in nto_cube_patterns):
            cube_source = os.path.join(orca_dir, file)
            cube_dest = os.path.join(trans_folder, file)
            move(cube_source, cube_dest)
            files_moved.append(file)
            print(f"   📦 Moved: {file} -> S{trans_num}/")
    
    if files_moved:
        print(f"   🗂️ Organised {len(files_moved)} files for S{trans_num}")
        return True
    else:
        print(f"   ⚠️  No files found to organize for S{trans_num}")
        return False



def parse_nto_for_selected_states(filepath, selected_excitations):
    """ Parse NTO data for selected excitation states and find minimum number of transitions needed to reach cumulative occupation >= 0.8. 
    Returns dict of selected states, along with corresponding properties
    e.g. {NTO 1: {"transition": int, "wavelength": float, "energy_ev": float, "fosc": float, "nto_transitions": list of dicts [{"from": int, "to": int, "occupation": float}, ...], "num_nto_needed": int, "total_occupation": float}, ...}. """
    
    with open(filepath, "r") as file:
        lines = file.readlines()
    
    all_transitions = []                                                                                    # Collect all selected excitation numbers.
    for peaks in selected_excitations.values():
        all_transitions.extend([p["transition"] for p in peaks])
    all_transitions = sorted(set(all_transitions))
    
    nto_results = {}
    
    for trans_num in all_transitions:                                                                       # For each selected excitation, parse NTO table for associated data.
        occupation_lines = []
        
        for i, line in enumerate(lines):
            if "NATURAL TRANSITION ORBITALS FOR STATE" in line and str(trans_num) in line:
                found_energy = False
                for j in range(i + 1, min(i + 50, len(lines))):
                    current_line = lines[j]
                    
                    if "NATURAL TRANSITION ORBITALS FOR STATE" in current_line:
                        break
                    
                    if not found_energy and current_line.strip().startswith("E="):
                        found_energy = True
                        continue
                    
                    if found_energy and "->" in current_line and ": n=" in current_line:
                        match = search(r"(\d+)a\s*->\s*(\d+)a\s*:\s*n=\s*([\d.]+)", current_line)           # Use regex to parse the occupation line.
                        if match:                            
                            occupation_lines.append({
                                "from": int(match.group(1)),                                                # Save "electron" orbital.
                                "to": int(match.group(2)),                                                  # Save "hole" orbital.
                                "occupation": float(match.group(3))                                         # Save "occupation" value, captures contribution of each NTO transition to overall excitation.
                            })
                                
                    if occupation_lines and current_line.strip() == "":
                            break        
                break
                    
        if not occupation_lines:
            print(f"Warning: NTO section not found for state {trans_num}")
            continue
        
        occupation_lines.sort(key=lambda x: x["occupation"], reverse=True)                                  # Sort by occupation (highest first).
        
        cumulative_sum = 0.0                                                                                # Find minimum number of transitions to reach sum >= 0.8.
        selected_nto = []
        for occ in occupation_lines:
            if cumulative_sum < 0.8:
                selected_nto.append(occ)
                cumulative_sum += occ["occupation"]
            else:
                break
        
        orig_transition = None                                                                              # Find the corresponding excitation from original data.
        for peaks in selected_excitations.values():
            for peak in peaks:
                if peak["transition"] == trans_num:
                    orig_transition = peak
                    break
            if orig_transition:
                break
        
        if orig_transition:
            nto_results[trans_num] = {"transition": trans_num, "wavelength": orig_transition["wavelength"], "energy_ev": orig_transition["energy_ev"], "fosc": orig_transition["fosc"], "nto_transitions": selected_nto, "num_nto_needed": len(selected_nto), "total_occupation": cumulative_sum}
    
    return nto_results



def parse_orca_tddft_absorption(filepath):
    """ Parse ORCA TD-DFT output file to extract absorption spectrum data.
    Return list of dicts [{"transition": int, "wavelength": float, "energy_ev": float, "fosc": float}, ...]. """
    
    with open(filepath, "r") as file:
        lines = file.readlines()

    start_line = None                                                                                               # Search for the absorption spectrum table, establish correct starting line for data.
    for i, line in enumerate(lines):
        if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line:
            for j in range(i + 1, min(i + 20, len(lines))):
                if "---" in lines[j] and j > i + 1:
                    start_line = j + 1
                    break
            break

    if start_line is None:
        raise ValueError("Absorption spectrum table not found")
    
    transitions = []
    for line in lines[start_line:]:
        if line.strip() == "" or "---" in line:
            break
        


        parts = line.strip().split()
        if len(parts) >= 8:
            try:                                                                                                    # The transition is in the format "0-1A", "->", "1-1A", "2.456315 (Energy (eV))", "19811.5 (Energy (cm-1))", "504.8 (Wavelength (nm))", "0.001791970 (fosc(D2))", ...
                trans_num = int(parts[2].split("-")[0])                                                             # Extract the transition number from third element (e.g., "1(a)" from "1(a)-1(b)A").
                wavelength = float(parts[5])                                                                        # Save sixth element (e.g., "504.8") as wavelength.
                energy_ev = photon_constant / wavelength                                                            # Calculate energy in eV.
                fosc = float(parts[6])                                                                              # Save seventh element (e.g., "0.001791970") as oscillation strength.
                
                transitions.append({"transition": trans_num, "wavelength": wavelength, "energy_ev": energy_ev, "fosc": fosc})
                    
            except (ValueError, IndexError):
                
                match = search(r"\d+-\S+\s*->\s*(\d+)-\S+\s+([\d.]+)\s+[\d.]+\s+([\d.]+)\s+([\d.]+)", line)         # Look for the pattern: number -> number, then 5 numbers.
                if match:
                    transitions.append({"transition": int(match.group(1)), "wavelength": float(match.group(3)), "energy_ev": photon_constant / float(match.group(3)), "fosc": float(match.group(4))})
    
    if not transitions:
        raise ValueError("No valid transition data found in the absorption spectrum table")
    
    return transitions



def run_multiwfn_generation(multiwfn_dir, orca_dir, filename, transition_numbers, nto_data):
    """ Run Multiwfn to generate NTOs and CDD files for multiple transitions. Processes one transition at a time to avoid file naming conflicts. """
    
    print(f"\n\n{'='*75}")
    print("📊 Step 5: Generating NTO folders for sorting")
    print(f"{'='*75}")
    
    workspace_info = setup_workspace(filename, transition_numbers)
    results = {}
    
    for i, trans_num in enumerate(transition_numbers, 1):
        success = generate_transition_files(multiwfn_dir, orca_dir, filename, trans_num, workspace_info, nto_data, i, len(transition_numbers))      # Generate NTO and CDD files.
        results[trans_num] = success
        if trans_num != transition_numbers[-1]:
            time.sleep(2)
    
    print(f"\n{'='*75}")
    print("🎉 Final Summary:")
    print(f"{'='*75}\n")
    

    for trans_num, success in results.items():
        print(f"   {'✅' if success else '❌'} Transition {trans_num}")
    
    print(f"\n   🗂️ Files organised in: {workspace_info['calc_folder']}")
    
    return results 



def save_selected_excitations(selected, nto_data, filename, calc_folder):
    """ Save selected excitations and NTO info to a text file in the Orbitals folder. """
    
    output_file = os.path.join(calc_folder, f"{filename}_selected_excitations.txt")
    
    with open(output_file, "w") as f:
        f.write("Selected Excitations from ORCA TD-DFT\n" + "=" * 60 + "\n\n")
        
        for category, peaks in selected.items():                                            # Print peak category (high energy or local).
            peaks_sorted = sorted(peaks, key=lambda x: x["transition"])
            f.write(f"{category.upper()}: {len(peaks_sorted)} peaks\n" + "-" * 40 + "\n")
            
            for peak in peaks_sorted:                                                       # Print excitation info.
                f.write(f"Transition {peak['transition']:3d}: {peak['wavelength']:6.1f} nm, {peak['energy_ev']:.3f} eV, fosc = {peak['fosc']:.3f}\n")
                
                if peak["transition"] in nto_data:                                          # Add NTO info if available.
                    nto_info = nto_data[peak["transition"]]
                    f.write(f"         NTOs needed: {nto_info['num_nto_needed']} (total occ: {nto_info['total_occupation']:.3f})\n")
                    for occ in nto_info["nto_transitions"]:
                        f.write(f"            {occ['from']}a -> {occ['to']}a : n = {occ['occupation']:.4f}\n")
                f.write("\n")
            f.write("\n")
    
    print(f"\n   💾 Results saved to '{output_file}'")



def select_excitation_peaks_by_criteria(transitions, en_thresh):
    """ Select excitations from absorption spectrum based on specific criteria:
    1. High intensity peaks (fosc >= 0.1)
    2. Local peaks (0.05 <= fosc < 0.1 and > 0.5eV away from any high intensity and local peaks)
    Returns dict {"high_intensity": high_intensity, "local_peaks": local_peaks}. """
    
    if not transitions:
        return {"high_intensity": [], "local_peaks": []}
    
    sorted_transitions = sorted(transitions, key=lambda x: x["wavelength"])
    
    
    high_intensity = [t for t in sorted_transitions if t["fosc"] >= 0.1]            # Select high energy peaks based on fosc strength.
    high_energies = [t["energy_ev"] for t in high_intensity]                        # Get energies of high intensity peaks for distance checking, need for next step.
    
    print(f"\n   🗂️ Selected {len(high_intensity)} high intensity peaks")
    
    
    local_candidates = []
    for t in sorted_transitions:
        if 0.05 <= t["fosc"] < 0.1:                                                 # Select local peaks based on fosc strength.
            if all(abs(t["energy_ev"] - he) > en_thresh for he in high_energies):   # Check energy gap between this local peak and other high energy peaks.
                local_candidates.append(t)                                          # Save local peak if all conditions met.
            else:
                print(f"   🗑️ Excluding local peak S{t['transition']}: within {en_thresh} eV of high intensity peak")
    
    local_peaks = []                                                                # Filter local candidates against each other (keep higher fosc when too close).
    if local_candidates:
        local_candidates.sort(key=lambda x: x["fosc"], reverse=True)
        
        for candidate in local_candidates:
            is_far = True
            for selected in local_peaks:
                if all(abs(candidate["energy_ev"] - selected["energy_ev"]) < en_thresh for selected in local_peaks):
                    is_far = False
                    print(f"   🗑️ Excluding local peak S{candidate['transition']} (fosc={candidate['fosc']:.4f}): "
                          f"within {en_thresh} eV of S{selected['transition']} (fosc={selected['fosc']:.4f})")
                    break
            if is_far:
                local_peaks.append(candidate)

    print(f"\n   🗂️ Selected {len(local_peaks)} local peaks")

    return {"high_intensity": high_intensity, "local_peaks": local_peaks}



def setup_workspace(filename, transition_numbers):
    """ Create organised workspace for NTO and CDD files. """
    
    current_dir = os.getcwd()                                                   # Create main Orbitals folder in the Code Directory if it doesn't exist.
    orbitals_base = os.path.join(current_dir, "Orbitals")
    os.makedirs(orbitals_base, exist_ok=True)
    
    calc_folder = os.path.join(orbitals_base, filename)                         # Create subfolder for this specific calculation.
    os.makedirs(calc_folder, exist_ok=True)
    
    transition_folders = {}                                                     # Create transition-specific folders within the calculation subfolder.
    for trans_num in transition_numbers:
        trans_folder = os.path.join(calc_folder, f"S{trans_num}")
        os.makedirs(trans_folder, exist_ok=True)
        transition_folders[trans_num] = trans_folder
    
    print("\n   🗂️ Workspace structure created:\n")
    print(f"    Parent Folder:        {orbitals_base}")
    print(f"    Calculation:          {calc_folder}")
    for trans_num, folder in transition_folders.items():
        if len(str(trans_num)) == 1:
            print(f"    S{trans_num}:                   {folder}")
        elif len(str(trans_num)) == 2:
            print(f"    S{trans_num}:                  {folder}")
    
    return {"orbitals_base": orbitals_base, "calc_folder": calc_folder, "transition_folders": transition_folders}



@contextmanager
def working_directory(path):
    """Context manager for temporarily changing directory"""
    original_dir = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(original_dir)
        
        
        