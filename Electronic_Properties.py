
# -*- coding: utf-8 -*-
"""
@author: dmalcolm
"""

import re

from Config import d_block_elements, p_block_elements




constant1 = 627.5       # Hartree to kcal/mol conversion factor.

# ============ HELPER FUNCTIONS ============

def check_frequency_file(content):
    """ Helper to check if file is a pure frequency calculation. Skips certain analysis if so. """
    
    return "GIBBS FREE ENERGY" in content and "OPT" not in content



def find_last_section(lines, section_name):
    """ Find the last occurrence of a section header in the file. Allows for access of the converged data specifically. """
    
    for i in range(len(lines) - 1, -1, -1):
        if section_name in lines[i]:
            return i
    return -1



def sort_dict_ntop(input_dict):
    """ Sorts a dictionary by values in ascending order (most negative to most positive). """

    return dict(sorted(input_dict.items(), key=lambda item: item[1]))



def sort_dict_pton(input_dict):
    """ Sorts a dictionary by values in descending order (most positive to most negative). """
    
    return dict(sorted(input_dict.items(), key=lambda item: item[1], reverse=True))



# ============ CHARGE EXTRACTION FUNCTIONS ============

def extract_hir_charges(orca_output_file, atom_labels):
    """ Extracts Hirshfeld charge values from ORCA output file.
    Returns dict {atom: charge}. """
    
    hirshfeld_charges = {}
    
    with open(orca_output_file, "r") as f:                                      # Early exit for frequency-only calculations.
        content = f.read()
        if check_frequency_file(content):
            return {}
    
    if not atom_labels:
        return hirshfeld_charges
    
    with open(orca_output_file, "r") as f:
        lines = f.readlines()
    
    
    last_section_pos = find_last_section(lines, "HIRSHFELD ANALYSIS")           # Find the last HIRSHFELD ANALYSIS section.
    if last_section_pos == -1:
        return hirshfeld_charges
    
    atom_counter = 0
    capture_data = False
    header_passed = False
    
    for line in lines[last_section_pos:]:
        if "HIRSHFELD ANALYSIS" in line:
            capture_data = True
            continue
        
        if capture_data and not header_passed:
            if "ATOM" in line and "CHARGE" in line and "SPIN" in line:
                header_passed = True
                continue
        
        if capture_data and header_passed:
            if "TOTAL" in line or (not line.strip() and atom_counter >= len(atom_labels)):
                break
            
            if line.strip() and line.strip()[0].isdigit():
                parts = line.split()
                try:
                    charge = float(parts[2])                                    # Charge is column index 2.
                    if atom_counter < len(atom_labels):
                        hirshfeld_charges[atom_labels[atom_counter]] = charge
                        atom_counter += 1
                    else:
                        break
                except (ValueError, IndexError) as e:
                    continue
    
    return hirshfeld_charges



def extract_mol_charge(filename):
    """ Looks for pattern: "* xyz m n", extracts "m" as the molecular charge  """
    
    try:
        with open(filename, "r") as f:
            for line in f:
                match = re.search(r"\*\s+xyz\s+([-+]?\d+)\s+\d+", line, re.IGNORECASE)      # Look for pattern: "* xyz -4 5" (charge is first number after xyz).
                if match:
                    charge = int(match.group(1))
                    return charge
    except Exception as e:
        print(f"Error reading molecular charge from {filename}: {e}")
    
    return None



def extract_qa_charges(orca_output_file, atom_labels):
    """ Extract Mulliken gross atomic charges (QA) from ORCA output.
    Returns dict {atom: charge}. """
    
    qa_charges = {}
    
    with open(orca_output_file, "r") as f:                                          # Early exit for frequency-only calculations.
        content = f.read()
        if check_frequency_file(content):
            return {}
    
    if not atom_labels:
        return qa_charges
    
    with open(orca_output_file, "r") as f:
        lines = f.readlines()
    
    last_section_pos = find_last_section(lines, "MAYER POPULATION ANALYSIS")        # Find the last MAYER POPULATION ANALYSIS section.
    if last_section_pos == -1:
        return qa_charges
    
    atom_counter = 0
    capture_data = False
    header_passed = False
    
    for line in lines[last_section_pos:]:
        if "MAYER POPULATION ANALYSIS" in line:
            capture_data = True
            continue
        
        if capture_data and not header_passed:
            if "ATOM" in line and "NA" in line and "QA" in line:
                header_passed = True
                continue
        
        if capture_data and header_passed:
            if not line.strip():
                continue
            
            if line.strip()[0].isdigit():                                           # Lines with data start with a digit.
                parts = line.split()
                try:
                    qa = float(parts[4])                                            # QA is column index 4 (0-based: ATOM, NA, ..., QA).
                    if atom_counter < len(atom_labels):
                        qa_charges[atom_labels[atom_counter]] = qa
                        atom_counter += 1
                    else:
                        break
                except (ValueError, IndexError) as e:
                    continue
    
    return qa_charges



# ============ SPIN FUNCTIONS ============

def extract_hir_pop(orca_output_file, atom_labels):
    """ Extracts Hirshfeld spin populations from ORCA output. 
    Returns dict {atom: spin}. """
    
    hirshfeld_spin = {}
    with open(orca_output_file, "r") as f:                                      # Early exit for frequency-only calculations.
        content = f.read()
        if check_frequency_file(content):
            return {}
        
    if not atom_labels:
        return hirshfeld_spin
    
    with open(orca_output_file, "r") as f:
        lines = f.readlines()
    
    last_section_pos = find_last_section(lines, "HIRSHFELD ANALYSIS")           # Find the last HIRSHFELD ANALYSIS section
    if last_section_pos == -1:
        return hirshfeld_spin
    
    atom_counter = 0
    capture_data = False
    header_passed = False
    
    for line in lines[last_section_pos:]:
        if "HIRSHFELD ANALYSIS" in line:
            capture_data = True
            continue
        
        if capture_data and not header_passed:
            if "ATOM" in line and "CHARGE" in line and "SPIN" in line:
                header_passed = True
                continue
        
        if capture_data and header_passed:
            if "TOTAL" in line or (not line.strip() and atom_counter >= len(atom_labels)):
                break
            
            if line.strip() and line.strip()[0].isdigit():
                parts = line.split()
                try:
                    spin = float(parts[3])                                      # Spin is column index 3
                    if atom_counter < len(atom_labels):
                        hirshfeld_spin[atom_labels[atom_counter]] = spin
                        atom_counter += 1
                    else:
                        break
                except (ValueError, IndexError) as e:
                    continue

    return hirshfeld_spin



def extract_loew_pop(orca_output_file, atom_labels):
    """ Extracts Loewdin spin populations from ORCA output. 
    Returns dict {atom: spin}. """
    
    loewdin_spin = {}
    
    with open(orca_output_file, "r") as f:                                      # Early exit for frequency-only calculations.
        content = f.read()
        if check_frequency_file(content):
            return {}
        
    if not atom_labels:
        return loewdin_spin
    
    with open(orca_output_file, "r") as f:
        lines = f.readlines()

    atom_counter = 0
    capture_data = False
    
    for line in lines:
        if "LOEWDIN ATOMIC CHARGES AND SPIN POPULATIONS" in line:
            capture_data = True
            continue
        
        if capture_data:
            if "-------------------------------------------" in line:
                continue
            
            if not line.strip():
                break
            
            if line.strip():
                cleaned_line = line.replace(":", "")                            # Remove all colons.
                parts = cleaned_line.split()
                try:
                    if len(parts) >= 4:
                        spin = float(parts[3])                                  # Spin is column index 3
                        if atom_counter < len(atom_labels):
                            loewdin_spin[atom_labels[atom_counter]] = spin
                            atom_counter += 1
                        else:
                            break
                    else:
                         continue
                except (ValueError, IndexError) as e:
                    continue
    
    return loewdin_spin



def extract_orb_spin_contribution(orca_output_file, atom_labels):
    """ Extract orbital-specific spin contributions for all atoms from ORCA output. 
    Returns dict {atom_label: {"orbitals": {orbital: spin, ...}, "total_orbital_spin": float}}. """
        
    orbital_spin_data = {}
        
    with open(orca_output_file, "r") as f:
        content = f.read()
        if check_frequency_file(content):
            return {}
        
    if not atom_labels:
        return orbital_spin_data
        
    with open(orca_output_file, "r") as f:
        lines = f.readlines()

    atom_to_index = {}
    element_map = {}
    for label in atom_labels:                                                                       # Extract the number from label like "Fe51" -> 51
        number = "".join(c for c in label if c.isdigit())
        element = "".join(c for c in label if c.isalpha())
        if number:
            atom_to_index[label] = int(number)
            element_map[int(number)] = element
        
    
    index_to_atom = {v: k for k, v in atom_to_index.items()}                                        # Create reverse mapping: index -> label
    
    d_orbital_order = ["dz2", "dxz", "dyz", "dx2y2", "dxy"]                                         # Define the orbital order for consistent display.
    p_orbital_order = ["px", "py", "pz"]
        
    in_section = False
    current_atom = None
    orbital_dict = {}
    
    for i, line in enumerate(lines):
        if "LOEWDIN REDUCED ORBITAL CHARGES AND SPIN POPULATIONS" in line:
            in_section = True
            continue
        
        if not in_section:
            continue
        
        if line.strip().startswith("-----") or line.strip().startswith("==="):                      # Look for atom lines like: "   51 Fes       :     0.013564"
            continue
        
        parts = line.strip().split()
        if parts and parts[0].isdigit() and len(parts) >= 2:

            if current_atom is not None and orbital_dict:
                
                label = index_to_atom.get(current_atom)                                             # Determine which orbital type to report.
                if label:
                    element = element_map.get(current_atom, "")
    
                    if element in d_block_elements:
                        target_orbitals = d_orbital_order
                    elif element in p_block_elements:
                        target_orbitals = p_orbital_order
                    else:
                        target_orbitals = ["s"]
                                
                    filtered_orbitals = {}                                                          # Find the maximum spin contribution among target orbitals.
                    total_spin = 0.0
                    for orb, spin in orbital_dict.items():
                        orb_clean = orb.replace("(", "").replace(")", "")
                        if orb_clean in target_orbitals or orb in target_orbitals:
                            filtered_orbitals[orb] = round(spin, 5)
                            total_spin += spin
                            
                    if filtered_orbitals:
                        orbital_spin_data[label] = {"orbitals": filtered_orbitals, "total_orbital_spin": round(total_spin, 5)}
                        
            current_atom = int(parts[0])                                                            # Start a new atom.
            
            try:
                colon_pos = line.find(":")
                if colon_pos != -1:
                    after_colon = line[colon_pos+1:].strip()
                    s_spin_parts = after_colon.split()                                              # The first number after the colon is the s-orbital spin.
                    if s_spin_parts:
                        s_spin = float(s_spin_parts[0])
                        orbital_dict["s"] = s_spin
            except (ValueError, IndexError):
                pass
            
            continue
                    
        if current_atom is not None and ":" in line and not line.strip().startswith("-----"):       # Collect orbital data for current atom, look for orbital lines: "      dz2     :     0.826355".
            parts = line.strip().split(":", 1)
            if len(parts) == 2:
                orbital_name = parts[0].strip()
                remainder = parts[1].strip()
                match = re.search(r"([-+]?\d*\.\d+)", remainder)
                if match:
                    spin_value = float(match.group(1))
                    orbital_dict[orbital_name] = spin_value
                    
        if not line.strip() and current_atom is not None:                                           # Stop when we hit a blank line after the section, save the last atoms' orbital data.
            if orbital_dict:
                label = index_to_atom.get(current_atom)
                if label:
                    element = element_map.get(current_atom, "")
                                
                    if element in d_block_elements:
                        target_orbitals = d_orbital_order
                    elif element in p_block_elements:
                        target_orbitals = p_orbital_order
                    else:
                        target_orbitals = ["s"]
                                
                    filtered_orbitals = {}
                    total_spin = 0.0
                    for orb, spin in orbital_dict.items():
                        orb_clean = orb.replace("(", "").replace(")", "")
                        if orb_clean in target_orbitals or orb in target_orbitals:
                            filtered_orbitals[orb] = round(spin, 5)
                            total_spin += spin
                            
                    if filtered_orbitals:
                        orbital_spin_data[label] = {
                            "orbitals": filtered_orbitals,
                            "total_orbital_spin": round(total_spin, 5)
                        }
            current_atom = None
            orbital_dict = {}
        
    return orbital_spin_data



# ============ ORBITAL EXTRACT FUNCTIONS ============

def extract_homo_lumo(file_path):
    """ Extracts HOMO/LUMO energies (eV) from non-spin-polarized ORCA output. 
    Returns dict {"homo": {"energy": float, "orbital": str}, "lumo": {"energy": float, "orbital": str}}. """
    
    results = {"homo": {"energy": None, "orbital": None}, "lumo": {"energy": None, "orbital": None}}

    try:
        with open(file_path, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        return results
    except Exception as e:
        return results

    orbital_starts = [i for i, line in enumerate(lines) if "ORBITAL ENERGIES" in line]          # Find the start of the last ORBITAL ENERGIES section.
    if not orbital_starts:
        return results
    
    last_orbital_start = orbital_starts[-1]
    
    orbitals = []
    for line in lines[last_orbital_start:]:
        if "ORBITAL ENERGIES" in line:
            continue
        if "Only the first 10 virtual orbitals" in line or "MOLECULAR ORBITALS" in line:
            break
        
        parsed = parse_orbital_line(line)
        if parsed:
            orbitals.append((parsed[1], parsed[2], parsed[0], parsed[3]))  # occ, energy, num, sym
    
    if orbitals:
        last_occupied = None
        first_unoccupied = None
        
        for occ, energy, orbital_num, sym in orbitals:
            if occ > 0:
                last_occupied = (energy, f"{orbital_num}{sym}")
            elif occ == 0 and first_unoccupied is None:
                first_unoccupied = (energy, f"{orbital_num}{sym}")
        
        if last_occupied:
            results["homo"]["energy"] = last_occupied[0]
            results["homo"]["orbital"] = last_occupied[1]
        if first_unoccupied:
            results["lumo"]["energy"] = first_unoccupied[0]
            results["lumo"]["orbital"] = first_unoccupied[1]
    
    return results



def extract_homo_lumo_from_spin_file(file_path):
    """ Extracts HOMO/LUMO energies (eV) from spin-polarized ORCA output.
    Spin up orbitals have default symmetry "α", spin down have default symmetry "β".
    Returns dict {"spin-α_homo": {"energy": float, "orbital": str}, "spin-α_lumo": {"energy": float, "orbital": str}, "spin-β_homo": {"energy": float, "orbital": str}, "spin-β_lumo": {"energy": float, "orbital": str}}. """

    
    results = {"spin-α_homo": {"energy": None, "orbital": None}, "spin-α_lumo": {"energy": None, "orbital": None}, "spin-β_homo": {"energy": None, "orbital": None}, "spin-β_lumo": {"energy": None, "orbital": None}}

    try:
        with open(file_path, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        return results
    except Exception as e:
        return results

    
    orbital_starts = [i for i, line in enumerate(lines) if "ORBITAL ENERGIES" in line]              # Find the start of the last ORBITAL ENERGIES section.
    if not orbital_starts:
        return results
    
    last_orbital_start = orbital_starts[-1]
    
    current_spin = None
    up_orbitals = []
    down_orbitals = []

    for line in lines[last_orbital_start:]:
        if "SPIN UP ORBITALS" in line:
            current_spin = "up"
            continue
        elif "SPIN DOWN ORBITALS" in line:
            current_spin = "down"
            continue
        elif "Only the first 10 virtual orbitals" in line or "MOLECULAR ORBITALS" in line:
            break
        
        if current_spin and line.strip():
            parsed = parse_orbital_line(line)
            if parsed:
                orbital_num, occ, energy, sym = parsed
                if current_spin == "up":
                    up_orbitals.append((occ, energy, orbital_num, sym))
                else:
                    down_orbitals.append((occ, energy, orbital_num, sym))
    
    if up_orbitals:                                                                                 # Find HOMO/LUMO for spin up.
        for occ, energy, num, sym in up_orbitals:
            if occ > 0:
                results["spin-α_homo"] = {"energy": energy, "orbital": f"{num}{sym}"}
            elif occ == 0 and results["spin-α_homo"]["energy"] is not None:
                results["spin-α_lumo"] = {"energy": energy, "orbital": f"{num}{sym}"}
                break
    
    if down_orbitals:                                                                               # Find HOMO/LUMO for spin down.
        for occ, energy, num, sym in down_orbitals:
            if occ > 0:
                results["spin-β_homo"] = {"energy": energy, "orbital": f"{num}{sym}"}
            elif occ == 0 and results["spin-β_homo"]["energy"] is not None:
                results["spin-β_lumo"] = {"energy": energy, "orbital": f"{num}{sym}"}
                break
    
    return results



def extract_orbitals_around_gap(file_path):
    """ Extracts HOMO/LUMO energies and 10 orbitals above LUMO and below HOMO from a non-spin-polarized file. 
    Returns dict {"homo": float, "lumo": float, "below_homo": {1: float, 2: float, ..., 10: float}, "above_lumo": {1: float, 2: float, ..., 10: float}}. """
    
    results = {"homo": None, "lumo": None, "below_homo": {}, "above_lumo": {}}

    try:
        with open(file_path, "r") as f:
            content = f.read()
    except Exception as e:
        return results

    marker = "ORBITAL ENERGIES"
    table_start = content.rfind(marker)          # last occurrence — works for both 1-table and 2-table outputs
    if table_start == -1:

        return results
    
    table_content = content[table_start:]
    found_header = False
    orbitals = []
    
    for line in table_content.split("\n"):
        if "NO   OCC" in line and "E(eV)" in line:
            found_header = True
            continue
        
        if not found_header or not line.strip():
            continue
        
        parsed = parse_orbital_line(line)
        if parsed:
            orbitals.append({
                "number": parsed[0], "occupation": parsed[1], "energy_ev": parsed[2]})
    
    homo_idx = None
    lumo_idx = None
    for i, orb in enumerate(orbitals):
        if orb["occupation"] == 2.0:
            results["homo"] = orb["energy_ev"]
            homo_idx = i
        elif orb["occupation"] == 0.0 and results["homo"] is not None and results["lumo"] is None:
            results["lumo"] = orb["energy_ev"]
            lumo_idx = i
    
    if results["homo"] is None or results["lumo"] is None:
        return results
    
    if homo_idx is not None:                                                            # Get 10 orbitals below HOMO.
        start = max(0, homo_idx - 10)
        for i, orb in enumerate(reversed(orbitals[start:homo_idx]), start=1):
            if i > 10:
                break
            results["below_homo"][i] = orb["energy_ev"]
    
    if lumo_idx is not None:                                                            # Get 10 orbitals above LUMO.
        end = min(len(orbitals), lumo_idx + 11)
        for i, orb in enumerate(orbitals[lumo_idx+1:end], start=1):
            if i > 10:
                break
            results["above_lumo"][i] = orb["energy_ev"]
    
    return results



def extract_orbitals_around_gap_from_spin_file(file_path):
    """ Extracts HOMO/LUMO from surrounding orbitals from spin-polarized file. 
    Returns dict {"spin-α_homo": float, "spin-α_lumo": float, "spin-β_homo": float, "spin-β_lumo": float, "spin-α_below_homo": {1: float, ...}, "spin-α_above_lumo": {1: float, ...}, "spin-β_below_homo": {1: float, ...}, "spin-β_above_lumo": {1: float, ...}}. """
    
    results = {"spin-α_homo": None, "spin-α_lumo": None, "spin-β_homo": None, "spin-β_lumo": None, "spin-α_below_homo": {}, "spin-α_above_lumo": {}, "spin-β_below_homo": {}, "spin-β_above_lumo": {}}

    try:
        with open(file_path, "r") as f:
            content = f.read()
    except Exception as e:
        return results


    orbital_start = content.rfind("ORBITAL ENERGIES")                                       # Find the ORBITAL ENERGIES section.
    if orbital_start == -1:
        return results
    
    next_section_positions = []
    for marker in ["MOLECULAR ORBITALS", "CARTESIAN COORDINATES", "TOTAL ENERGY"]:          # Find the end of this section, look for common section markers after the orbital table or file end.
        pos = content.find(marker, orbital_start + 1)
        if pos != -1:
            next_section_positions.append(pos)
    
    next_section = min(next_section_positions) if next_section_positions else len(content)
    orbital_block = content[orbital_start:next_section]
    
    up_start = orbital_block.find("SPIN UP ORBITALS")
    down_start = orbital_block.find("SPIN DOWN ORBITALS")
    
    if up_start != -1 and down_start != -1 and down_start > up_start:
        up_block = orbital_block[up_start:down_start]                                       # Extract spin up block (from start of spin up to start of spin down tables).
        down_block = orbital_block[down_start:]                                             # Extract spin down block (from start of spin down table till end).
        
        up_data = process_orbital_block(up_block)
        down_data = process_orbital_block(down_block)                   
        
        results["spin-α_homo"] = up_data["homo"]
        results["spin-α_lumo"] = up_data["lumo"]
        results["spin-α_below_homo"] = up_data["below_homo"]
        results["spin-α_above_lumo"] = up_data["above_lumo"]
        
        results["spin-β_homo"] = down_data["homo"]
        results["spin-β_lumo"] = down_data["lumo"]
        results["spin-β_below_homo"] = down_data["below_homo"]
        results["spin-β_above_lumo"] = down_data["above_lumo"]

    return results



def parse_orbital_line(line):
    """ Parse a single orbital energy line from ORCA output. 
    Returns list [orbital number, occupation, energy (eV), symmetry (spin up or down)]. """
    
    parts = line.split()
    if len(parts) < 4 or not parts[0].isdigit():
        return None
    
    try:
        orb_num = int(parts[0])
        occ = float(parts[1])
        energy_ev = float(parts[3])
        symmetry = parts[4] if len(parts) >= 5 else "a"

        return (orb_num, occ, energy_ev, symmetry)
    except (ValueError, IndexError):
        return None



def process_orbital_block(block):
    """ Process a block of orbital data (either spin up or down). Extracts frontier and surrounding orbitals, as well as their energy values. """
    
    orbitals = []
    homo_index = None
    lumo_index = None
    results = {"homo": None, "lumo": None, "below_homo": {}, "above_lumo": {}}
    
    for line in block.split("\n"):
        if (not line.strip() or "NO" in line or "OCC" in line or 
            "Only the first" in line or "---" in line):
            continue
        
        parsed = parse_orbital_line(line)
        if parsed:
            orbitals.append({"number": parsed[0], "occupation": parsed[1], "energy_ev": parsed[2], "symmetry": parsed[3]})
    
    for i, orb in enumerate(orbitals):
        if orb["occupation"] > 0.5:                 # Occupied orbital.
            results["homo"] = orb["energy_ev"]
            homo_index = i
        elif orb["occupation"] < 0.5 and results["homo"] is not None and results["lumo"] is None:
            results["lumo"] = orb["energy_ev"]
            lumo_index = i
            break
    
    if homo_index is not None:                      # Get 10 orbitals below HOMO.
        start = max(0, homo_index - 10)
        for i, orb in enumerate(reversed(orbitals[start:homo_index]), start=1):
            if i > 10:
                break
            results["below_homo"][i] = orb["energy_ev"]
    
    if lumo_index is not None:                      # Get 10 orbitals above LUMO.
        end = min(len(orbitals), lumo_index + 11)
        for i, orb in enumerate(orbitals[lumo_index+1:end], start=1):
            if i > 10:
                break
            results["above_lumo"][i] = orb["energy_ev"]
    
    return results



# ============ ORBITAL JOIN FUNCTIONS ============

def extract_homo_lumo_join(file_path):
    """ Extracts HOMO/LUMO energies, auto-detects spin polarization. """
    
    spin_results = extract_homo_lumo_from_spin_file(file_path)                                                                         # Establish if spin-polarisation present.
    
    spin_data_found = (spin_results["spin-α_homo"]["energy"] is not None or spin_results["spin-β_homo"]["energy"] is not None)         # Confirm if spin-polarized data was found (check energy values).
    
    if spin_data_found:
        spin_results["spin_polarized"] = True
        return spin_results
    
    
    else:                                                                                                                              # If no spin data found, report back as non-polarized.
        non_spin_results = extract_homo_lumo(file_path)
        non_spin_results["spin_polarized"] = False
        return non_spin_results



def extended_extract_homo_lumo_join(file_path):
    """ Extract HOMO-10 to LUMO+10, auto-detecting spin polarization. """
    
    spin_results = extract_orbitals_around_gap_from_spin_file(file_path)                                                               # Establish if spin-polarisation present.
    
    # If no spin data found, try non-polarized
    if (spin_results["spin-α_homo"] is not None or spin_results["spin-β_homo"] is not None):
        return spin_results
    else:
        return extract_orbitals_around_gap(file_path)



# ============ ELECTRONIC ENERGY FUNCTIONS ============

def extract_energy_values(filename):
    """ Extract electronic energy and thermal corrections (if present) from ORCA output.
    Returns dict {{"calculation_type": str, "electronic_energy_eh": float, "thermal_correction_eh": float, "electronic_energy_kcal": float, "thermal_correction_kcal": float}}. """
    
    energy_data = {"calculation_type": None, "electronic_energy_eh": None, "thermal_correction_eh": None, "electronic_energy_kcal": None, "thermal_correction_kcal": None}
    
    try:
        with open(filename, "r") as f:
            lines = f.readlines()
            content = "".join(lines)
    except Exception as e:
        return {}

    for i, line in enumerate(lines):
        if "INPUT FILE" in line and i + 3 < len(lines):
            keyword_line = lines[i+3]
            if re.search(r"\bFREQ\b", keyword_line, re.IGNORECASE):
                energy_data["calculation_type"] = "Frequency"
            elif re.search(r"\bNUMFREQ\b", keyword_line, re.IGNORECASE):
                energy_data["calculation_type"] = "Numerical Frequency"
            elif re.search(r"\bSP\b", keyword_line, re.IGNORECASE):
                energy_data["calculation_type"] = "Single Point"
            elif re.search(r"\bOPT\b", keyword_line, re.IGNORECASE):
                energy_data["calculation_type"] = "Optimization"
            break
    
    if energy_data["calculation_type"] in ["Frequency", "Numerical Frequency"]:                                     # Determine calculation type.
        
        gibbs_match = re.search(r"Final Gibbs free energy\s+\.+\s+([-\d\.]+)\s+Eh", content)                        # Extract Gibbs free energy.
        if gibbs_match:
            energy_data["electronic_energy_eh"] = float(gibbs_match.group(1))
            energy_data["electronic_energy_kcal"] = energy_data["electronic_energy_eh"] * constant1
        
        thermal_match = re.search(r"G-E\(el\)\s+\.+\s+([-\d\.]+)\s+Eh\s+([-\d\.]+)\s+kcal/mol", content)            # Extract thermal correction (in kcal/mol).
        if thermal_match:
            energy_data["thermal_correction_eh"] = float(thermal_match.group(1))
            energy_data["thermal_correction_kcal"] = float(thermal_match.group(2))
    
    elif energy_data["calculation_type"] in ["Single Point", "Optimization"]:
        
        sp_energies = re.findall(r"FINAL SINGLE POINT ENERGY\s+([-\d\.]+)", content)                                # Find all single point energies and take the last one.
        
        if sp_energies:
            energy_data["electronic_energy_eh"] = float(sp_energies[-1])  # Last occurrence
            energy_data["electronic_energy_kcal"] = energy_data["electronic_energy_eh"] * constant1
    
    else:
        return {}

    return energy_data



def extract_vib_freq(filename, calculation_type):
    """ Extracts vibrational frequencies from ORCA frequency calculation.
    Returns dict {{"all_freq": {mode: frequency, ...}, "imag_freq": {mode: negative_frequency, ...}, "has_imag": bool, "num_imag": int}}. """
    
    
    if calculation_type not in ["Frequency", "Numerical Frequency"]:                                                                   # Only extract frequencies for frequency calculations
        return {}    
    
    try:
        with open(filename, "r") as f:
            content = f.read()
    except Exception as e:
        return {}
    
    freq_section_match = re.search(r"-+\s*VIBRATIONAL FREQUENCIES\s*-+.*?(?=-+\s*\w+\s*-+)", content, re.DOTALL | re.IGNORECASE)       # Find the VIBRATIONAL FREQUENCIES section
    
    if not freq_section_match:
        return {}
    
    freq_section = freq_section_match.group(0)
    
    freq_lines = re.findall( r"\s+(\d+):\s+([-\d\.]+)\s+cm\*\*-1", freq_section)                                                       # Extract all frequency lines (format: "    164:    3715.48 cm**-1")
    
    if not freq_lines:
        return {}
    
    all_freq = {}
    imag_freq = {}

    for mode_num, freq_value in freq_lines:
        freq = float(freq_value)
        mode = int(mode_num)
        all_freq[mode] = freq
        if freq < 0:
            imag_freq[mode] = freq
            
    freq_results = {"all_freq": all_freq, "imag_freq": imag_freq, "has_imag": len(imag_freq) > 0, "num_imag": len(imag_freq)}
    
    if freq_results["has_imag"]:
        print(f"""WARNING: Found {freq_results["num_imag"]} imaginary frequency(ies)!""")
        for mode, freq in imag_freq.items():
            print(f"  Mode {mode}: {freq:.2f} cm**-1")
    
    return freq_results


