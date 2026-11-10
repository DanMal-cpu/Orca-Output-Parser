# -*- coding: utf-8 -*-
"""
@author: dmalcolm
"""

from math import sqrt
from re import finditer




# ============ HELPER FUNCTIONS ============

def atom_sort_key(atom_name):
    """ Generate a sort key for atom labels based on element and number 
    e.g., "O20" -> ("O", 20). """
    
    element = "".join(filter(str.isalpha, atom_name))
    number = "".join(filter(str.isdigit, atom_name))
    return (element, int(number) if number else 0)



def bond_key(bond_string):
    """ Generate a sort key for bond strings.
    e.g., "Fe1-O2" -> (("Fe", 1), ("O", 2)). """
    
    if "-" in bond_string:
        a1, a2 = bond_string.split("-")
        return (atom_sort_key(a1), atom_sort_key(a2))
    return (("", 0), ("", 0))



def clean_atom(atom_str):
    """ Standardize atom format, including removing parentheses and fixing "0O" -> "O0" atom format. """
    
    clean = atom_str.replace("(", "").replace(")", "").replace(" ", "")
    parts = clean.split("-")
    
    if len(parts) == 2:
        if parts[0].isdigit() and parts[1].isalpha():
            return f"{parts[1]}{parts[0]}"                  # "0O" → "O0"
        elif parts[1].isdigit() and parts[0].isalpha():
            return f"{parts[0]}{parts[1]}"                  # "Mo40" → "Mo40"
    return clean



def get_element(atom):
    """ Extract element symbol from atom label, handling both 1 and 2 letter element symbols. 
    e.g. "O" from "O20". """
    
    return atom[:2] if len(atom) > 1 and atom[1].islower() else atom[0]



def sort_bond_dict(bond_dict):
    """ Sort bond dictionary for consistent ordering. """
    
    return dict(sorted(bond_dict.items(), key=lambda x: bond_key(x[0])))



def standardize_bond(atom1, atom2):
    """ Standardize bond order: oxygen always second, otherwise alphabetical.
    e.g. "H1-O2", "Fe3-H1", ... """
    
    o1 = atom1.startswith("O")
    o2 = atom2.startswith("O")
    
    if o1 and not o2:
        return f"{atom2}-{atom1}"                               # Swap to put oxygen second.
    elif not o1 and o2:
        return f"{atom1}-{atom2}"                               # Oxygen already second.
    else:
        return f"{min(atom1, atom2)}-{max(atom1, atom2)}"       # Alphabetical.



# ============ ORCA SPECIFIC FUNCTIONS ============

def calculate_bond_lengths(coordinates, bl_cutoff):
    """ Calculate bond lengths between all atom pairs, filtering by element-specific cutoffs.
    Returns dict {bond name: float, ...}. """
    
    n_atoms = len(coordinates)
    bond_lengths = {}
    
    pair_cutoffs = {("O", "O"): 1.5}                                # Special cutoffs for specific element pairs.
    
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            a1, x1, y1, z1 = coordinates[i]
            a2, x2, y2, z2 = coordinates[j]
            
            bond = standardize_bond(a1, a2)                         # Standardize bond designation (oxygen always second).
            
            dist = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)       # Calculate distance.
            dist = round(dist, 4)
            
            e1, e2 = get_element(a1), get_element(a2)               # Get element-specific cutoff.
            pair = tuple(sorted((e1, e2)))
            
            if pair in pair_cutoffs:
                max_allowed = pair_cutoffs[pair]
            else:
                max_allowed = min(bl_cutoff.get(e1, 3.0), bl_cutoff.get(e2, 3.0))
            
            if dist <= max_allowed:                                 # Store if within cutoff.
                bond_lengths[bond] = dist
    
    return sort_bond_dict(bond_lengths)



def extract_and_label_coordinates(filename):
    """ Extract Cartesian coordinates and label atoms sequentially.
    Return list of tuples [(atom_id, x, y, z), ...]. """
    
    with open(filename, "r") as f:
        lines = f.readlines()
    
    header_indices = [i for i, line in enumerate(lines) if "CARTESIAN" in line.upper() and "ANGSTROEM" in line.upper()]                # Find last CARTESIAN COORDINATES header.
    
    if not header_indices:
        return []
    
    last_header = header_indices[-1]
    
    coordinates = []
    atom_counter = 0
    
    for line in lines[last_header + 1:]:
        line = line.strip()
        
        if line.startswith("---") or "COORDINATES" in line.upper():
            continue
        
        if coordinates and not line:
            break
        
        parts = line.split()
        if len(parts) >= 3:
            coords = []
            for part in reversed(parts):
                try:
                    coords.insert(0, float(part))
                    if len(coords) == 3:
                        break
                except ValueError:
                    pass
            
            if len(coords) == 3:
                atom_name = next((p for p in parts if not p.replace(".", "").replace("-", "").isdigit()), "X")
                labeled_atom = f"{atom_name}{atom_counter}"
                coordinates.append((labeled_atom, coords[0], coords[1], coords[2]))
                atom_counter += 1
    
    return coordinates



def extract_bond_orders(filename):
    """ Extract Mayer bond orders from ORCA output file.
    Returns dict {bond name: float, ...}. """
    
    try:
        with open(filename, "r") as f:
            content = f.read()
    except Exception as e:
        print(f"Error reading file: {e}")
        return {}
    
    start_idx = content.find("Mayer bond orders larger than")
    if start_idx == -1:
        print("Error: Could not find bond orders section.")
        return {}
    
    pattern = r"B\(([^,]+),\s*([^)]+)\)\s*:\s*([\d.]+)"
    bond_orders = {}
    
    for match in finditer(pattern, content[start_idx:]):
        atom1, atom2, value = match.groups()
        atom1, atom2 = clean_atom(atom1), clean_atom(atom2)
        bond = standardize_bond(atom1, atom2)
        bond_orders[bond] = float(value)
    
    return sort_bond_dict(bond_orders) if bond_orders else {}



def store_atom_coordinates(atom_data):
    """ Store coordinates in separate dictionaries for easy access. 
    Return list of tuples ([atom_id, ...], {atom_id: x, ...}, {atomd_id: y, ...}, {atom_id: z, ...}). """
    
    atom_ids = []
    x, y, z = {}, {}, {}
    
    for atom_id, xc, yc, zc in atom_data:
        atom_ids.append(atom_id)
        x[atom_id] = xc
        y[atom_id] = yc
        z[atom_id] = zc
    
    return atom_ids, x, y, z
