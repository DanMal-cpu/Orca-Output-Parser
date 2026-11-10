# -*- coding: utf-8 -*-
"""
@author: dmalcolm
"""

import numpy as np

from collections import defaultdict
from re import IGNORECASE, search

from Config import class_groups
from Molecule_Classification import class_order, class_to_param, format_superscript_charge, generate_molecule_formula




class OrcaParameterStorage:
    """ Storage and formatting class for ORCA output parameters. This class manages all extracted parameters and provides methods for formatting them for display or export. """
    
    def __init__(self):
        """ Initialize empty parameter storage. """
        
        self.parameters = {}                                    
        self.requested_parameters = set()                       
 
    
    def get_calc_type(self):
        """ Retrieve calculation type """
        calc_type = "Unknown"
        energy_data = self.parameters.get("energy_data", {})
        if "calculation_type" in energy_data:
            calc_type = energy_data["calculation_type"]
        
        return calc_type
    
        
    # ============ FORMATTING METHODS ============
    
    def format_atom_class(self, atom_class_data):
        """ Format atom classification data in a hierarchical, organized format. """
        
        if not atom_class_data:
            return "No atom classification data available.\n\n"
        
        total_atoms = len(atom_class_data)
         
        class_to_atoms = {}                                                                     # Reverse mapping: classification -> list of atom labels.
        for atom, classification in atom_class_data.items():
            if classification not in class_to_atoms:
                class_to_atoms[classification] = []
            class_to_atoms[classification].append(atom)
        
        def natural_atom_key(atom):
            """ Sort atoms within each classification numerically (e.g. O3, O19, O39 - not O19, O3, O39). """
            atom_str = str(atom)
            digits = "".join(c for c in atom_str if c.isdigit())
            return (int(digits) if digits else 0)

        for cls in class_to_atoms:
            class_to_atoms[cls].sort(key=natural_atom_key)
         
        groups_to_display = []
        displayed_classes = set()
        
        for group_name, group_info in class_groups.items():
            classifications = group_info["classifications"]
            subgroup_data = []
            
            for cls in classifications:
                if cls in class_to_atoms:
                    atoms = class_to_atoms[cls]
                    subgroup_data.append((cls, atoms))
                    displayed_classes.add(cls)
            
            if not subgroup_data:
                continue
            
            # Inclusion logic:
            if "OVERALL" in group_name:                                                         # Either a) always include groups with "OVERALL" in the name.
                include = True
            else:                                                                               # b) include groups that are NOT subsets of any OVERALL group.
                is_subset_of_overall = False
                for other_name, other_info in class_groups.items():
                    if "OVERALL" in other_name and other_name != group_name:                    # Check if all our classifications are in that OVERALL group.
                        other_classifications = set(other_info["classifications"])
                        if all(cls in other_classifications for cls in classifications):
                            is_subset_of_overall = True
                            break
                include = not is_subset_of_overall
            
            if include:
                display_name = group_name.replace(" - OVERALL", "")
                groups_to_display.append((display_name, subgroup_data))
        
        for classification, atoms in class_to_atoms.items():                                    # Add any classifications that weren't in any group.
            if classification not in displayed_classes:
                groups_to_display.append((classification, [(classification, atoms)]))
        

        class_width = 27
        output = []
        
        output.append(f"""\n\n{"CLASS":<{class_width}} COUNT   ATOMS""")
        
        divider_idx = len(output)
        output.append("")                                                                       # Placeholder for divider.
        
        first = True
        for display_name, subgroup_data in groups_to_display:
            if not first:
                output.append("")
            first = False
            
            total_count = sum(len(atoms) for _, atoms in subgroup_data)
            
            if len(subgroup_data) > 1:
                output.append(f"{display_name:<{class_width}} {total_count:>5}   ")             # Group with subgroups.
                for cls, atoms in subgroup_data:
                    atom_string = ", ".join(atoms)
                    cls_display = self.format_subclass_name(cls, display_name)
                    full_name = f"  - {cls_display}"
                    output.append(f"{full_name:<{class_width}} {len(atoms):>5}   {atom_string}")
            else:
                cls, atoms = subgroup_data[0]                                                   # Single classification.
                atom_string = ", ".join(atoms)
                output.append(f"{display_name:<{class_width}} {total_count:>5}   {atom_string}")
        
        # Adjust the divider to match the longest line
        max_len = 0                                                                             # Find the longest line (excluding header and footer and the divider itself).
        for i, line in enumerate(output):
            if i != divider_idx and not line.startswith("=") and line.strip():                  # Skip header/footer/divider.
                max_len = max(max_len, len(line))
        
        divider_len = max(max_len + 1, 50)
        output[divider_idx] = "-" * divider_len
        
        return "\n".join(output)


    
    def format_bond_table(self, bond_lengths, bond_orders):
        """ Format bond lengths and bond orders into a single table, ordered element first.
        Returns formatted table with columns: Bond Name, Bond Length, Bond Order. """
        
        if not bond_lengths or not bond_orders:
            return "No bond data available"
        
        common_bonds = set(bond_lengths.keys()) & set(bond_orders.keys())               # Get all unique bonds that appear in bond_lengths or bond_orders dictionary.
        if not common_bonds:
            return "No bond data available"
        
        def extract_numbers(bond_str):
            """Extract numeric parts from bond string for sorting"""
            bond_str = str(bond_str)
            if "-" not in bond_str:
                return (0, 0)
    
            parts = bond_str.split("-")
            if len(parts) < 2:
                return (0, 0)
            
            def to_int(part):
                digits = "".join(c for c in part if c.isdigit())
                return int(digits) if digits else 0
            
            return (to_int(parts[0]), to_int(parts[1]))
        
        element_groups = defaultdict(list)                                              # Group bonds by first element.
        for bond in common_bonds:
            bond_str = str(bond)
            element = bond_str.split("-")[0] if "-" in bond_str else "other"
            element_groups[element].append(bond)
        
        lines = []
        lines.append(f"""{"BOND NAME":<12} {"BOND LENGTH":<12} {"BOND ORDER":<12}""")
        lines.append("-" * 45)
        
        sorted_elements = sorted(element_groups.keys())
        
        for i, element in enumerate(sorted_elements):
            bonds = element_groups[element]
            if not bonds:
                continue
        
            sorted_bonds = sorted(bonds, key=extract_numbers)
            for bond in sorted_bonds:
                length = bond_lengths.get(bond, "N/A")
                order = bond_orders.get(bond, "N/A")
                
                if isinstance(length, float) and isinstance(order, float):
                    lines.append(f"{bond:<12} {length:<12.4f} {order:<12.4f}")
                else:
                    lines.append(f"{bond:<12} {str(length):<12} {str(order):<12}")
            
            if i < len(sorted_elements) - 1:
                lines.append("")
        
        return "\n".join(lines)

 

    def format_class_data(self, classification_data):
        """ Format classification data as a simple atom-to-class list, with atoms grouped by classification and sorted numerically.
        Returns formatted string listing each atom and its classification. """

        if not classification_data:
            return "No data available"
        
        class_groups = defaultdict(list)                                                                            # Group atoms by their classification.
        for atom, atom_class in classification_data.items():
            class_groups[atom_class].append(atom)
        
        formatted_lines = []
        
        def add_class_group(atom_class, atoms):
            """Helper to format and add a classification group"""
            if not atoms:
                return
                
            sorted_atoms = sorted(atoms, key=lambda x: int("".join(c for c in str(x) if c.isdigit())))              # Sort atoms numerically.
            
            max_atom_len = max(len(str(atom)) for atom in sorted_atoms)                                             # Find maximum atom ID length for alignment.
            
            for atom in sorted_atoms:                                                                               # Add each atom with alignment.
                spacing = max_atom_len - len(str(atom)) + 2
                formatted_lines.append(f"""{atom}:{" " * spacing}{atom_class}""")
            
            formatted_lines.append("")                                                                              # Add spacing between groups.
        
        for atom_class in class_order:                                                                              # Process predefined classes in order.
            if atom_class in class_groups:
                add_class_group(atom_class, class_groups[atom_class])
                del class_groups[atom_class]                                                                        # Remove from dict so it won't appear in "remaining" section.
        
        for atom_class in sorted(class_groups.keys()):                                                              # Handle any remaining classes.
            add_class_group(atom_class, class_groups[atom_class])
        
        return "\n".join(formatted_lines)    
    
    
    
    def format_coord(self, coords):
        """ Format atomic coordinates in a readable table. """
        
        if not coords or not all(isinstance(coord, tuple) and len(coord) == 4 for coord in coords):
            return str(coords)
        
        max_atom_id_len = max(len(str(coord[0])) for coord in coords)               # Find the maximum atom ID length for proper alignment.
        
        return "\n".join(f"{atom_id:<{max_atom_id_len}}   {x:9.6f}  {y:9.6f}  {z:9.6f}" for atom_id, x, y, z in coords)
    
    
    
    def format_dict(self, data):
        """ Format dictionary data with aligned values """
        
        if not data:
            return "No data available"
        
        element_groups = defaultdict(list)                  # Group by first element for organization.
        for key, value in data.items():
            key_str = str(key)
            first_element = key_str.split("-")[0] if "-" in key_str else "other"
            element_groups[first_element].append((key, value))
                
        formatted_lines = []
        sorted_elements = sorted(element_groups.keys())
        
        for element in sorted_elements:
            bonds = element_groups[element]
                
            def sort_bond_key(bond_item):
                """ Sort bonds numerically within each element group. """
                bond, value = bond_item
                bond_str = str(bond)
                if "-" in str(bond):
                    return (0, 0)
                    
                parts = bond_str.split("-")
                if len(parts) < 2:                          # Handle case where split doesn't give 2 parts.
                    return (0, 0)
            
                def get_number(part):
                    digits = "".join(c for c in part if c.isdigit())
                    return int(digits) if digits else 0
            
                return (get_number(parts[0]), get_number(parts[1]))
            
            sorted_bonds = sorted(bonds, key=sort_bond_key)
            max_bond_len = max(len(str(bond)) for bond, _ in sorted_bonds)
            
            for bond, value in sorted_bonds:                # Add each bond in this element group.
                bond_str = str(bond)
                if isinstance(value, float):
                    formatted_line = f"""{bond_str}:{" " * (max_bond_len - len(bond_str) + 2)}{value:.4f}"""
                else:
                    formatted_line = f"{bond_str}: {value}"
                formatted_lines.append(formatted_line)
            
            if element != sorted_elements[-1]:              # Add spacing between different element groups.
                formatted_lines.append("")
        
        return "\n".join(formatted_lines)    
    
  
    
    def format_elec_energy_data(self, energy_data):
        """ Format electronic energy data for display. Includes appropriate units. """
        
        lines = []
    
        if not energy_data:
            return "No energy data available"
        
        def safe_format(value, precision):
            """Format Eh values with 8 decimal places"""
            if value is None:
                return "N/A"
            try:
                if isinstance(value, (int, float)):
                    return f"{value:.{precision}f}"
                return str(value)
            except (ValueError, TypeError):
                return str(value)


        ee_eh = safe_format(energy_data.get("electronic_energy_eh"), 8)                     # Electronic Energy
        ee_kcal = safe_format(energy_data.get("electronic_energy_kcal"), 2)
        lines.append(f"Electronic Energy:     {ee_eh} Eh")
        lines.append(f"Electronic Energy:     {ee_kcal} kcal/mol")
        lines.append("")
        
        tc_eh = safe_format(energy_data.get("thermal_correction_eh"), 8)                    # Thermal Correction
        tc_kcal = safe_format(energy_data.get("thermal_correction_kcal"), 2)
        lines.append(f"Thermal Correction:    {tc_eh} Eh")
        lines.append(f"Thermal Correction:    {tc_kcal} kcal/mol")
        
        if "enthalpy_eh" in energy_data:                                                    # Optional: Add enthalpy if available
            h_eh = safe_format(energy_data.get("enthalpy_eh"), 8)
            h_kcal = safe_format(energy_data.get("enthalpy_kcal"), 2)
            lines.append("")
            lines.append(f"Enthalpy (H):          {h_eh} Eh")
            lines.append(f"Enthalpy (H):          {h_kcal} kcal/mol")
        
        if "gibbs_eh" in energy_data:                                                       # Optional: Add Gibbs free energy if available
            g_eh = safe_format(energy_data.get("gibbs_eh"), 8)
            g_kcal = safe_format(energy_data.get("gibbs_kcal"), 2)
            lines.append("")
            lines.append(f"Gibbs Free Energy (G): {g_eh} Eh")
            lines.append(f"Gibbs Free Energy (G): {g_kcal} kcal/mol")
        
        if "calculation_type" in energy_data and energy_data["calculation_type"]:           # Add calculation type if available.
            lines.append("")
            lines.append(f"""Calculation Type:      {energy_data["calculation_type"]}""")
        
        return "\n".join(lines)
 
 

    def format_frontier_data(self, frontier_data, extended_data):
        """ Format frontier orbital data for display. Handles both spin-polarized and non-polarized cases. """

        lines = []
        
        if frontier_data:
            lines.append("FRONTIER ORBITALS:")
            lines.append("-" * 40)
            
            spin_polarized = frontier_data.get("spin_polarized", False)
            
            
            if spin_polarized:
                alpha_homo = frontier_data.get("spin-α_homo", {})                                                       # Handle spin-polarized case (α and β spins).
                alpha_lumo = frontier_data.get("spin-α_lumo", {})
                beta_homo = frontier_data.get("spin-β_homo", {})
                beta_lumo = frontier_data.get("spin-β_lumo", {})
                
                alpha_homo_energy = alpha_homo.get("energy", "N/A")                                                     # Alpha (up) spin.
                alpha_homo_orbital = alpha_homo.get("orbital", "N/A")
                alpha_lumo_energy = alpha_lumo.get("energy", "N/A")
                alpha_lumo_orbital = alpha_lumo.get("orbital", "N/A")
                
                beta_homo_energy = beta_homo.get("energy", "N/A")                                                       # Beta (down) spin.
                beta_homo_orbital = beta_homo.get("orbital", "N/A")
                beta_lumo_energy = beta_lumo.get("energy", "N/A")
                beta_lumo_orbital = beta_lumo.get("orbital", "N/A")
                
                
                if isinstance(alpha_homo_energy, float) and isinstance(alpha_lumo_energy, float):                       # Calculate alpha spin gap.
                    alpha_gap = alpha_lumo_energy - alpha_homo_energy
                else:
                    alpha_gap = "N/A"
                    
                if isinstance(beta_homo_energy, float) and isinstance(beta_lumo_energy, float):                         # Calculate beta spin gap.
                    beta_gap = beta_lumo_energy - beta_homo_energy
                else:
                    beta_gap = "N/A"
                
                
                lines.append("ALPHA SPIN:")                                                                             # Format alpha spin data.
                if isinstance(alpha_homo_energy, float):
                    lines.append(f"  HOMO(α):    {alpha_homo_energy:8.4f} eV  (orbital {alpha_homo_orbital})")
                else:
                    lines.append(f"  HOMO(α):    {alpha_homo_energy:>8} eV  (orbital {alpha_homo_orbital})")
                
                if isinstance(alpha_lumo_energy, float):
                    lines.append(f"  LUMO(α):    {alpha_lumo_energy:8.4f} eV  (orbital {alpha_lumo_orbital})")
                else:
                    lines.append(f"  LUMO(α):    {alpha_lumo_energy:>8} eV  (orbital {alpha_lumo_orbital})")
                    
                if isinstance(alpha_gap, float):
                    lines.append(f"  Gap(α):     {alpha_gap:8.4f} eV")
                else:
                    lines.append(f"  Gap(α):     {alpha_gap:>8} eV")
                
                lines.append("")
                
                
                lines.append("BETA SPIN:")                                                                              # Format beta spin data.
                if isinstance(beta_homo_energy, float):
                    lines.append(f"  HOMO(β):    {beta_homo_energy:8.4f} eV  (orbital {beta_homo_orbital})")
                else:
                    lines.append(f"  HOMO(β):    {beta_homo_energy:>8} eV  (orbital {beta_homo_orbital})")
                
                if isinstance(beta_lumo_energy, float):
                    lines.append(f"  LUMO(β):    {beta_lumo_energy:8.4f} eV  (orbital {beta_lumo_orbital})")
                else:
                    lines.append(f"  LUMO(β):    {beta_lumo_energy:>8} eV  (orbital {beta_lumo_orbital})")
                    
                if isinstance(beta_gap, float):
                    lines.append(f"  Gap(β):     {beta_gap:8.4f} eV")
                else:
                    lines.append(f"  Gap(β):     {beta_gap:>8} eV")
                    
                lines.append("")
                    
                    
            else:                                                                                                       # Handle non-spin-polarized case.
                homo_energy = frontier_data.get("homo", {}).get("energy", "N/A")
                homo_orbital = frontier_data.get("homo", {}).get("orbital", "N/A")
                lumo_energy = frontier_data.get("lumo", {}).get("energy", "N/A") 
                lumo_orbital = frontier_data.get("lumo", {}).get("orbital", "N/A")
                
                if isinstance(homo_energy, float) and isinstance(lumo_energy, float):
                    gap = lumo_energy - homo_energy
                else:
                    gap = "N/A"
                
                
                if isinstance(homo_energy, float):                                                                      # Format HOMO line.
                    lines.append(f"HOMO:        {homo_energy:8.4f} eV  (orbital {homo_orbital})")
                else:
                    lines.append(f"HOMO:        {homo_energy:>8} eV  (orbital {homo_orbital})")
                
                if isinstance(lumo_energy, float):                                                                      # Format LUMO line.
                    lines.append(f"LUMO:        {lumo_energy:8.4f} eV  (orbital {lumo_orbital})")
                else:
                    lines.append(f"LUMO:        {lumo_energy:>8} eV  (orbital {lumo_orbital})")
                    
                if isinstance(gap, float):                                                                              # Format HOMO-LUMO gap line.
                    lines.append(f"Gap:         {gap:8.4f} eV")
                else:
                    lines.append(f"Gap:         {gap:>8} eV")
                    
                lines.append("")
            
            lines.append(f"""Spin polarized: {"Yes" if spin_polarized else "No"}""")
            lines.append("")
        
        
        if extended_data:                                                                                               # Extended frontier orbitals in two-column layout.
            lines.append("")
            lines.append("")

            
            has_alpha_below = "spin-α_below_homo" in extended_data                                                      # Check if extended data has spin-polarized structure.
            has_alpha_above = "spin-α_above_lumo" in extended_data  
            has_beta_below = "spin-β_below_homo" in extended_data
            has_beta_above = "spin-β_above_lumo" in extended_data
            
            extended_spin_polarized = (has_alpha_below or has_alpha_above or 
                                     has_beta_below or has_beta_above)
            
            if extended_spin_polarized:                                                                                 # Handle spin-polarized extended data with defaults.
                
                alpha_below_homo = extended_data.get("spin-α_below_homo", {})
                alpha_above_lumo = extended_data.get("spin-α_above_lumo", {})
                beta_below_homo = extended_data.get("spin-β_below_homo", {})
                beta_above_lumo = extended_data.get("spin-β_above_lumo", {})
                alpha_homo_ext = extended_data.get("spin-α_homo")
                alpha_lumo_ext = extended_data.get("spin-α_lumo")
                beta_homo_ext = extended_data.get("spin-β_homo")
                beta_lumo_ext = extended_data.get("spin-β_lumo")
                
                has_alpha_data = (alpha_below_homo or alpha_above_lumo or alpha_homo_ext is not None or alpha_lumo_ext is not None)
                has_beta_data = (beta_below_homo or beta_above_lumo or beta_homo_ext is not None or beta_lumo_ext is not None)
                
                if has_alpha_data:                                                                                      # Alpha spin extended orbitals
                    
                    lines.append("EXTENDED FRONTIER ORBITALS (ALPHA SPIN):")
                    lines.append("-" * 60)
                    lines.append(f"""{"ORBITAL":<12} {"ENERGY (eV)":<12} | {"ORBITAL":<12} {"ENERGY (eV)":<12}""")
                    lines.append("-" * 60)
                    
                    alpha_left_column = []
                    alpha_right_column = []
                    
                    for i in range(10, 0, -1):                                                                          # Alpha spin orbitals below HOMO.
                        energy = alpha_below_homo.get(i)
                        if energy is not None:
                            if isinstance(energy, float):
                                alpha_left_column.append((f"HOMO-{i}(α)", f"{energy:8.4f}"))
                            else:
                                alpha_left_column.append((f"HOMO-{i}(α)", f"{energy:>8}"))
                    
                    if alpha_homo_ext is not None:                                                                      # Alpha HOMO and LUMO
                        if isinstance(alpha_homo_ext, float):
                            alpha_left_column.append(("HOMO(α)", f"{alpha_homo_ext:8.4f}"))
                        else:
                            alpha_left_column.append(("HOMO(α)", f"{alpha_homo_ext:>8}"))
                    if alpha_lumo_ext is not None:
                        if isinstance(alpha_lumo_ext, float):
                            alpha_left_column.append(("LUMO(α)", f"{alpha_lumo_ext:8.4f}"))
                        else:
                            alpha_left_column.append(("LUMO(α)", f"{alpha_lumo_ext:>8}"))
                    
                    for i in range(1, 11):                                                                              # Alpha spin orbitals above LUMO.
                        energy = alpha_above_lumo.get(i)
                        if energy is not None:
                            if isinstance(energy, float):
                                alpha_right_column.append((f"LUMO+{i}(α)", f"{energy:8.4f}"))
                            else:
                                alpha_right_column.append((f"LUMO+{i}(α)", f"{energy:>8}"))
                    
                    max_alpha_lines = max(len(alpha_left_column), len(alpha_right_column))
                    for i in range(max_alpha_lines):
                        left_orbital = ""
                        left_energy = ""
                        right_orbital = ""
                        right_energy = ""
                        
                        if i < len(alpha_left_column):
                            left_orbital, left_energy = alpha_left_column[i]
                        
                        if i < len(alpha_right_column):
                            right_orbital, right_energy = alpha_right_column[i]
                        
                        line = f"{left_orbital:<12} {left_energy:<12} | {right_orbital:<12} {right_energy:<12}"
                        lines.append(line.rstrip())
                    
                    lines.append("")
                
                
                if has_beta_data:                                                                                       # Beta spin extended orbitals
                    
                    lines.append("")
                    lines.append("")
                    lines.append("EXTENDED FRONTIER ORBITALS (BETA SPIN):")
                    lines.append("-" * 60)
                    lines.append(f"""{"ORBITAL":<12} {"ENERGY (eV)":<12} | {"ORBITAL":<12} {"ENERGY (eV)":<12}""")
                    lines.append("-" * 60)
                    
                    beta_left_column = []
                    beta_right_column = []
                    
                    for i in range(10, 0, -1):                                                                          # Beta spin orbitals below HOMO.
                        energy = beta_below_homo.get(i)
                        if energy is not None:
                            if isinstance(energy, float):
                                beta_left_column.append((f"HOMO-{i}(β)", f"{energy:8.4f}"))
                            else:
                                beta_left_column.append((f"HOMO-{i}(β)", f"{energy:>8}"))
                    
                    if beta_homo_ext is not None:                                                                       # Beta HOMO and LUMO
                        if isinstance(beta_homo_ext, float):
                            beta_left_column.append(("HOMO(β)", f"{beta_homo_ext:8.4f}"))
                        else:
                            beta_left_column.append(("HOMO(β)", f"{beta_homo_ext:>8}"))
                    if beta_lumo_ext is not None:
                        if isinstance(beta_lumo_ext, float):
                            beta_left_column.append(("LUMO(β)", f"{beta_lumo_ext:8.4f}"))
                        else:
                            beta_left_column.append(("LUMO(β)", f"{beta_lumo_ext:>8}"))
                    
                    for i in range(1, 11):                                                                              # Beta spin orbitals above LUMO.
                        energy = beta_above_lumo.get(i)
                        if energy is not None:
                            if isinstance(energy, float):
                                beta_right_column.append((f"LUMO+{i}(β)", f"{energy:8.4f}"))
                            else:
                                beta_right_column.append((f"LUMO+{i}(β)", f"{energy:>8}"))
                    
                    max_beta_lines = max(len(beta_left_column), len(beta_right_column))
                    for i in range(max_beta_lines):
                        left_orbital = ""
                        left_energy = ""
                        right_orbital = ""
                        right_energy = ""
                        
                        if i < len(beta_left_column):
                            left_orbital, left_energy = beta_left_column[i]
                        
                        if i < len(beta_right_column):
                            right_orbital, right_energy = beta_right_column[i]
                        
                        line = f"{left_orbital:<12} {left_energy:<12} | {right_orbital:<12} {right_energy:<12}"
                        lines.append(line.rstrip())
                    
            else:
                lines.append(f"""{"Orbital":<12} {"Energy (eV)":<12} | {"Orbital":<12} {"Energy (eV)":<12}""")          # Handle non-spin-polarized extended data.
                lines.append("-" * 60)
                
                below_homo = extended_data.get("below_homo", {})
                above_lumo = extended_data.get("above_lumo", {})
                homo_ext = extended_data.get("homo")
                lumo_ext = extended_data.get("lumo")
                
                has_data = (below_homo or above_lumo or homo_ext is not None or lumo_ext is not None)                   # Check if we have any data to display.
                
                if has_data:

                    left_column = []
                    right_column = []
                    
                    for i in range(10, 0, -1):                                                                          # Below HOMO orbitals (HOMO-10 to HOMO-1) - left column.
                        energy = below_homo.get(i)
                        if energy is not None:
                            if isinstance(energy, float):
                                left_column.append((f"HOMO-{i}", f"{energy:8.4f}"))
                            else:
                                left_column.append((f"HOMO-{i}", f"{energy:>8}"))
                    
                    if homo_ext is not None:                                                                            # HOMO and LUMO - left column.
                        if isinstance(homo_ext, float):
                            left_column.append(("HOMO", f"{homo_ext:8.4f}"))
                        else:
                            left_column.append(("HOMO", f"{homo_ext:>8}"))
                    if lumo_ext is not None:
                        if isinstance(lumo_ext, float):
                            left_column.append(("LUMO", f"{lumo_ext:8.4f}"))
                        else:
                            left_column.append(("LUMO", f"{lumo_ext:>8}"))
                    
                    for i in range(1, 11):                                                                              # Above LUMO orbitals (LUMO+1 to LUMO+10) - right column.
                        energy = above_lumo.get(i)
                        if energy is not None:
                            if isinstance(energy, float):
                                right_column.append((f"LUMO+{i}", f"{energy:8.4f}"))
                            else:
                                right_column.append((f"LUMO+{i}", f"{energy:>8}"))
                    
                    max_lines = max(len(left_column), len(right_column))
                    
                    for i in range(max_lines):
                        left_orbital = ""
                        left_energy = ""
                        right_orbital = ""
                        right_energy = ""
                        
                        if i < len(left_column):
                            left_orbital, left_energy = left_column[i]
                        
                        if i < len(right_column):
                            right_orbital, right_energy = right_column[i]
                        
                        line = f"{left_orbital:<12} {left_energy:<12} | {right_orbital:<12} {right_energy:<12}"
                        lines.append(line.rstrip())
                else:
                    lines.append("No extended frontier orbital data available")
        
        return "\n".join(lines)

   
 
    def format_hir_table(self, hir_values, sorted_hir):
        """ Format Hirshfeld charge values using the generic two-column table. """
        
        return self.format_two_column_table(hir_values, sorted_hir, "HIR VALUE", "SORTED HIR VALUE")
    
  
    
    def format_loew_table(self, loew_values, sorted_loew):
        """ Format Löwdin spin populations using the generic two-column table. """
        
        return self.format_two_column_table(loew_values, sorted_loew, "LÖWDIN SPIN", "SORTED SPIN")
  
    
  
    def format_parameter(self, param_name, value):
        """ Format a parameter value for display based on its type and name. """
        
        if isinstance(value, dict):
            if param_name == "atom_class":
                return self.format_class_data(value)
            elif param_name == "energy_data":
                return self.format_elec_energy_data(value)
            elif param_name == "vib_freq":
                return self.format_vib_freq(value)
            elif "bond" in param_name:
                return self.format_dict(value) 
            elif "qa" in param_name:
                if "sort_qa" in param_name:
                    return self.format_qa_charges(value, is_sorted=True)
                else:
                    return self.format_qa_charges(value, is_sorted=False)
            elif "loewdin" in param_name:
                if "loewdin_sort" in param_name:
                    return self.format_loew_table(value, is_sorted=True)
                else:
                    return self.format_loew_table(value, is_sorted=False)
            else:
                return self.format_dict(value)
        
        elif isinstance(value, list):                                       # Deals with all list parameters.
            if param_name == "coords":
                return self.format_coord(value)
            else:
                return "\n".join(str(item) for item in value)
        else:
            return str(value)
    
    
  
    def format_qa_table(self, qa_charges, sorted_qa):
        """ Format QA (Mulliken) values using the generic two-column table. """
        
        return self.format_two_column_table(qa_charges, sorted_qa, "QA VALUE", "SORTED QA VALUE")
  
   
  
    def format_qa_charges(self, qa_data, is_sorted):
        """ Format QA (Mulliken) values. """
        
        if not qa_data:
            return "No data available"
        
        def extract_number(atom):
            return int("".join(c for c in str(atom) if c.isdigit()))
        
        def format_line(atom, value, max_len):
            atom_str = str(atom)
            padding = " " * (max_len - len(atom_str) + 2)
            value_str = f"{value:.4f}" if isinstance(value, float) else str(value)
            return f"{atom_str}:{padding}{value_str}"
        
        lines = []
        
        if is_sorted:                                                                           # Sort all atoms by their value.
            max_atom_len = max(len(str(atom)) for atom in qa_data)
            sorted_items = sorted(qa_data.items(), key=lambda x: x[1])
            
            for atom, value in sorted_items:
                lines.append(format_line(atom, value, max_atom_len))
        else:
            element_groups = defaultdict(list)
            for atom, value in qa_data.items():
                atom_str = str(atom)
                element = "".join([c for c in atom_str if c.isalpha()])
                element_groups[element].append((atom, value))
            
            for element in sorted(element_groups.keys()):
                atoms = element_groups[element]
                sorted_atoms = sorted(atoms, key=lambda x: extract_number(x[0]))                # Sort atoms numerically.
                max_atom_len = max(len(str(atom)) for atom, _ in sorted_atoms)
                
                for atom, value in sorted_atoms:
                    lines.append(format_line(atom, value, max_atom_len))
                
                lines.append("")                                                                # Add gap after each element group.
            
            if len(element_groups) == 1:                                                        # Calculate and add mean and standard deviation for single-element groups.
                values = [v for v in qa_data.values() if isinstance(v, (int, float))]
                if values:
                    lines.append(f"Mean: {np.mean(values):.4f}")
                    lines.append(f"Stdev: {np.std(values):.4f}")
        
        return "\n".join(lines)
 
    
 
    def format_subclass_name(self, cls, group_name):
        """ Format subclassification name based on the group pattern.
        
        Examples:
            - "Terminal Oxygen" with group "OXYGEN - TERMINAL (M=O)" -> "M=O Terminal"
            - "M-M Bridging Oxygen" with group "OXYGEN - BRIDGING (M-O-M)" -> "M-O-M Bridging"
            - "Regular Framework Metal" with group "FRAMEWORK METAL (M-O-M)" -> "M-O-M Regular"
            - "TM-Bridging Framework Metal" with group "FRAMEWORK METAL (M-O-TM)" -> "M-O-TM TM-Bridging" """
        
        pattern_match = search(r"\(([^)]+)\)", group_name)
        pattern = pattern_match.group(1) if pattern_match else None
        
        cls_display = cls
        for suffix in [" Framework Metal", " Oxygen", " Metal"]:                # Remove common suffixes.
            if cls_display.endswith(suffix):
                cls_display = cls_display.replace(suffix, "")
        
        if pattern:
            if "OXYGEN" in group_name:
                # For oxygen: pattern + " " + type
                # e.g., "M=O" + " Terminal" -> "M=O Terminal"
                return f"{pattern} {cls_display}"
            elif "FRAMEWORK METAL" in group_name:
                # For framework metals: pattern + " " + type
                # e.g., "M-O-M" + " Regular" -> "M-O-M Regular"
                return f"{pattern} {cls_display}"
        
        return cls_display                                                      # If no pattern or not matching, return cleaned classification name.
    


    def format_two_column_table(self, left_values, right_values, left_label, right_label):
        """ Generic function to format any two-column table. """
        
        if not left_values or not right_values:
            return f"No {left_label} data available"
            
        atom_width = 6
        value_width = 14
            
        def fmt_line(atom, value):
            atom_str = str(atom)
            atom_fmt = f"{atom_str}".ljust(atom_width)
                
            if isinstance(value, float):
                value_fmt = f"{value:>{value_width}.4f}"
            else:
                value_fmt = f"{value:>{value_width}}"
                
            return f"{atom_fmt}{value_fmt}"
            
        if isinstance(right_values, list):
            right_values = dict(right_values)
        elif not isinstance(right_values, dict):
            try:
                right_values = dict(right_values)
            except:
                return f"Error: right_values is type {type(right_values)}"
            
        left_lines = []
        elem_groups = defaultdict(list)
            
        for atom, val in left_values.items():
            elem = "".join(c for c in str(atom) if c.isalpha())
            elem_groups[elem].append((atom, val))
            
        for elem in sorted(elem_groups):
            atoms = sorted(elem_groups[elem], key=lambda x: int("".join(c for c in str(x[0]) if c.isdigit()) or 0))
            for atom, val in atoms:
                left_lines.append((atom, val, True))  # Has data
            left_lines.append((None, None, False))    # Empty line
            
        if left_lines and not left_lines[-1][2]:
            left_lines.pop()
            
        right_lines = [(a, v, True) for a, v in right_values.items()]
            
        sep = "     | "
        lines = []
        lines.append(f"""{"ATOM":<{atom_width}}{left_label:>{value_width}}{sep}{"ATOM":<{atom_width}}{right_label:>{value_width}}""")
        lines.append("-" * (2*(atom_width+value_width) + len(sep) + 3))
            
        left_idx = right_idx = 0
        while left_idx < len(left_lines) or right_idx < len(right_lines):
            if left_idx < len(left_lines):
                a1, v1, has_data = left_lines[left_idx]
                left = fmt_line(a1, v1) if has_data else " " * (atom_width + value_width)
                left_idx += 1
            else:
                left = " " * (atom_width + value_width)
                
            if right_idx < len(right_lines):
                a2, v2, _ = right_lines[right_idx]
                right = fmt_line(a2, v2)
                right_idx += 1
            else:
                right = " " * (atom_width + value_width)
                
            lines.append(f"{left}{sep}{right}")
            
        return "\n".join(lines)


 
    def format_vib_freq(self, vib_freq_data):
        """ Format vibrational frequency data for display. Scans for any imaginary frequencies """

        if not vib_freq_data:
            return "No vibrational frequency data available"  
   
        frequencies = []
        lines = []
        
        if isinstance(vib_freq_data, dict):
            if "all_freq" in vib_freq_data:
                freq_data = vib_freq_data["all_freq"]
                if isinstance(freq_data, dict):
                    frequencies = [freq_data[k] for k in sorted(freq_data.keys())]
                elif isinstance(freq_data, list):
                    frequencies = freq_data
                else:
                    return f"""Unexpected type for "all_freq": {type(freq_data)}"""
            else:
                return f"Frequency data not found. Available keys: {list(vib_freq_data.keys())}"  
        elif isinstance(vib_freq_data, list):
            frequencies = vib_freq_data
        else:
            return f"Unexpected data type: {type(vib_freq_data)}"
        
        if not frequencies:
            return "No frequency data available"
        
        num_of_freqs = len(frequencies)
        
        num_of_columns = 4
        rows_per_column = (num_of_freqs + num_of_columns - 1) // num_of_columns
        
        mode_width = 5
        freq_width = 12
        sep = " | "
        
        lines.append("=" * 90)
        lines.append(f"Vibrational Frequencies (Total: {num_of_freqs})")
        
        if isinstance(vib_freq_data, dict):                                             # Raise warning about imaginary frequencies if found.
            has_imag = vib_freq_data.get("has_imag", False)
            num_imag = vib_freq_data.get("num_imag", 0)
            if has_imag and num_imag > 0:
                lines.append(f"Warning: {num_imag} imaginary frequencies found")
        
        lines.append("=" * 90)
        lines.append("")
        

        col_headers = []
        for col in range(num_of_columns):
            col_headers.append(f"""{"Mode":<{mode_width}} {"Freq (cm⁻¹)":>{freq_width}}""")
        lines.append(sep.join(col_headers))
        lines.append("-" * (len(sep.join(col_headers))))
        
        col_data = []
        for col in range(num_of_columns):                                               # Organize frequencies into columns.
            start_idx = col * rows_per_column
            end_idx = min(start_idx + rows_per_column, num_of_freqs)
            if isinstance(frequencies, list):
                col_freqs = frequencies[start_idx:end_idx]
            else:
                freq_list = list(frequencies.values()) if isinstance(frequencies, dict) else list(frequencies)
                col_freqs = freq_list[start_idx:end_idx]
            col_data.append(col_freqs)
           
        max_rows = max(len(col) for col in col_data)                                    # Find the maximum number of rows needed. 

        for row in range(max_rows):                                                     # Create data rows, this will display ALL frequencies.
            row_parts = []
            for col in range(num_of_columns):
                if row < len(col_data[col]):
                    freq = col_data[col][row]
                    mode_num = (col * rows_per_column) + row + 1
                    
                    if isinstance(freq, (int, float)):
                        row_parts.append(f"{mode_num:<{mode_width}} {freq:>{freq_width}.2f}")
                    else:
                        row_parts.append(f"{mode_num:<{mode_width}} {str(freq):>{freq_width}}")
                else:
                    row_parts.append(f"""{"":<{mode_width}} {"":>{freq_width}}""")
            lines.append(sep.join(row_parts))
        
        lines.append("")                                                                # Add footer with statistics.
        lines.append("-" * 90)
        numeric_freqs = [f for f in frequencies if isinstance(f, (int, float))]
        if numeric_freqs:
            lines.append(f"Frequency Range: {min(numeric_freqs):.2f} - {max(numeric_freqs):.2f} cm⁻¹")
        
        return "\n".join(lines)   
 

    # ============ SUMMARY GENERATION METHODS ============ 

    def generate_basic_summary(self):
        """ Generate a formatted basic summary block for the trimmed file.
        Returns a str containing:
            - Compound formula (with charge)
            - Calculation type
            - Electronic energy (Eh and kcal/mol)
            - Spin polarization status
            - Frontier orbital gap (α and β if spin-polarized) """
        
        lines = []
        lines.append("")
        lines.append("===============")
        lines.append("BASIC SUMMARY")
        lines.append("===============")
        lines.append("")
    
        charge = self.parameters.get("charge")
        formula = generate_molecule_formula(self.parameters.get("atom_class"), charge)              # Molecular formula
        
        if formula:
            charge_str = format_superscript_charge(charge) if charge else ""
            lines.append(f"COMPOUND:          [{formula}]{charge_str}\n")
        else:
            lines.append("COMPOUND:          N/A\n")
    
        filename = self.parameters.get("filename", "")
        calc_type = "Unknown"
        energy_data = self.parameters.get("energy_data", {})
        if "calculation_type" in energy_data:
            calc_type = energy_data["calculation_type"]

        calc_match = search(r"(SP|Opt|Freq|NUMFREQ)", filename, IGNORECASE)
        calc_abbr = calc_match.group(0).upper() if calc_match else ""
        
        if calc_abbr:
            lines.append(f"CALCULATION:       {calc_type} ({calc_abbr})\n")
        else:
            lines.append(f"CALCULATION:       {calc_type}\n")

        ee_eh = energy_data.get("electronic_energy_eh")                                             # Electronic energy
        ee_kcal = energy_data.get("electronic_energy_kcal")
        if ee_eh is not None and ee_kcal is not None:
            lines.append(f"ELECTRONIC ENERGY: {ee_eh:.2f} Eh  |  {ee_kcal:,.2f} kcal/mol\n")
        else:
            lines.append("ELECTRONIC ENERGY: N/A\n")
    
        frontier = self.parameters.get("frontier", {})                                              # Spin polarization
        spin_polarized = frontier.get("spin_polarized", False)  
        lines.append(f"""SPIN POLARIZATION: {"Yes (α/β asymmetry detected)" if spin_polarized else "No"}\n""")
    
        if spin_polarized:                                                                          # Frontier gap
            alpha_homo = frontier.get("spin-α_homo", {}).get("energy")
            alpha_lumo = frontier.get("spin-α_lumo", {}).get("energy")
            beta_homo = frontier.get("spin-β_homo", {}).get("energy")
            beta_lumo = frontier.get("spin-β_lumo", {}).get("energy")
            if None not in (alpha_homo, alpha_lumo, beta_homo, beta_lumo):
                alpha_gap = alpha_lumo - alpha_homo
                beta_gap = beta_lumo - beta_homo
                lines.append(f"FRONTIER GAP:      α: {alpha_gap:.2f} eV  |  β: {beta_gap:.2f} eV\n")
            else:
                lines.append("FRONTIER GAP:      N/A\n")
        else:
            homo = frontier.get("homo", {}).get("energy")
            lumo = frontier.get("lumo", {}).get("energy")
            if homo is not None and lumo is not None:
                gap = lumo - homo
                lines.append(f"FRONTIER GAP:      {gap:.2f} eV\n")
            else:
                lines.append("FRONTIER GAP:      N/A\n")
    
        lines.append("")                                                                            # Blank line before next section
        return "\n".join(lines)   

 

    def generate_bond_summary(self):
        """
        Generate a bond summary table using the already-classified data.
        Returns a str containing:
            - Bond type
            - Count
            - Mean length
            - Mean order """
        
        atom_class = self.parameters.get("atom_class")
        if not atom_class:
            return "No atom classification data available"

        bond_params = [key for key in self.parameters.keys() if key.endswith("_bond_lengths")]
        
        if not bond_params:
            return "No bond data available"
        
        output = []
        output.append("=" * 35)
        output.append("BOND SUMMARY")
        output.append("=" * 35)
        output.append("")
        
        output.append(f"""{"BOND TYPE":<20} {"COUNT":>6} {"MEAN LENGTH":>12} {"MEAN ORDER":>11}""")
        output.append("-" * 53)
        
        bond_data = defaultdict(lambda: {"lengths": [], "orders": [], "count": 0})
        
        for param_name in sorted(bond_params):

            prefix = param_name.replace("_bond_lengths", "")
            
            base_type = self.get_bond_type_display(prefix)
            if base_type is None:
                continue
            
            bond_lengths = self.parameters.get(param_name, {})
            bond_orders_key = param_name.replace("_bond_lengths", "_bond_orders")
            bond_orders = self.parameters.get(bond_orders_key, {})
            
            if not bond_lengths:
                continue
            
            if prefix in ["mbridtm_O", "hetbridm_O", "hetbridtm_O", "tmtolig_O"]:
                for bond in bond_lengths.keys():
                    if bond not in bond_orders:
                        continue
                    
                    if "-" in bond:
                        atom1, atom2 = bond.split("-")
                        cls1 = atom_class.get(atom1, "")
                        cls2 = atom_class.get(atom2, "")
                        
                        if "Oxygen" in cls1:                                                                    # Determine which atom is oxygen.
                            other_cls = cls2
                            other_atom = atom2
                        elif "Oxygen" in cls2:
                            other_cls = cls1
                            other_atom = atom1
                        else:
                            continue
                        
                        if prefix == "mbridtm_O":

                            if "Framework Metal" in other_cls:
                                bond_type = "TM-(O-M)"
                            elif "TM-Substitute" in other_cls:
                                bond_type = "M-(O-TM)"
                            else:
                                continue
                                
                        elif prefix == "hetbridm_O":

                            if "Framework Metal" in other_cls:
                                bond_type = "X-(O-M)"
                            elif "Heteroatom" in other_cls:
                                bond_type = "M-(O-X)"
                            else:
                                continue
                                
                        elif prefix == "hetbridtm_O":

                            if "Heteroatom" in other_cls:
                                bond_type = "TM-(O-X)"
                            elif "TM-Substitute" in other_cls:
                                bond_type = "X-(O-TM)"
                            else:
                                continue

                        elif prefix == "tmtolig_O":

                            if "TM-Substitute" in other_cls:
                                bond_type = "Lig-(O-TM)"
                            elif "Ligand" in other_cls:
                                bond_type = "TM-(O-Lig)"
                            else:
                                continue

                        bond_data[bond_type]["lengths"].append(bond_lengths[bond])
                        bond_data[bond_type]["orders"].append(bond_orders[bond])
                        bond_data[bond_type]["count"] += 1
            
            else:
                for bond in bond_lengths.keys():
                    if bond not in bond_orders:
                        continue
                    bond_data[base_type]["lengths"].append(bond_lengths[bond])
                    bond_data[base_type]["orders"].append(bond_orders[bond])
                    bond_data[base_type]["count"] += 1
        
        if not bond_data:
            return "No bond data available"
                                                                                                                
        sorted_bond_types = sorted(bond_data.keys(), key=lambda x: (-bond_data[x]["count"], x))                 # Sort by count (descending)
        
        for bond_type in sorted_bond_types:
            data = bond_data[bond_type]
            count = data["count"]
            mean_length = np.mean(data["lengths"]) if data["lengths"] else 0
            mean_order = np.mean(data["orders"]) if data["orders"] else 0
            output.append(f"{bond_type:<20} {count:>6} {mean_length:>12.3f} {mean_order:>11.3f}")
        
        output.append("")
        
        return "\n".join(output)
    
        
    
    def generate_charge_summary(self, charge_type):
        """ Generate a charge summary table organized by classification groups.
        Returns formatted charge summary with count, mean, and standard deviation for each classification group. """
        
        if charge_type == "qa":
            charge_values = self.parameters.get("qa_charges")
            charge_label = "QA Mulliken Charges"
        #elif charge_type == "hir":
        #    charge_values = self.parameters.get("hir_charges")
        #    charge_label = "Hirshfeld Charges"        
        else:
            return f"Unknown charge type: {charge_type}"
        
        atom_class = self.parameters.get("atom_class")
        
        if not atom_class or not charge_values:
            return "No charge data available"
        
        output = []
        output.append("=" * 40)
        output.append(f"CHARGE SUMMARY ({charge_label})")
        output.append("=" * 40)
        output.append("")
        
        class_to_atoms = {}
        for atom, classification in atom_class.items():
            if classification not in class_to_atoms:
                class_to_atoms[classification] = []
            class_to_atoms[classification].append(atom)
        
        groups_to_display = []
        displayed_classes = set()
        
        for group_name, group_info in class_groups.items():
            classifications = group_info["classifications"]
            subgroup_data = []
            
            for cls in classifications:
                if cls in class_to_atoms:
                    atoms = class_to_atoms[cls]
                    values = []
                    for atom in atoms:
                        if atom in charge_values:
                            values.append(charge_values[atom])
                    subgroup_data.append((cls, atoms, values))
                    displayed_classes.add(cls)
            
            if not subgroup_data:
                continue
            
            if "OVERALL" in group_name:
                include = True
            else:
                is_subset_of_overall = False
                for other_name, other_info in class_groups.items():
                    if "OVERALL" in other_name and other_name != group_name:
                        other_classifications = set(other_info["classifications"])
                        if all(cls in other_classifications for cls in classifications):
                            is_subset_of_overall = True
                            break
                include = not is_subset_of_overall
            
            if include:
                display_name = group_name.replace(" - OVERALL", "")
                groups_to_display.append((display_name, subgroup_data))
        
        for classification, atoms in class_to_atoms.items():
            if classification not in displayed_classes:
                values = [charge_values[atom] for atom in atoms if atom in charge_values]
                groups_to_display.append((classification, [(classification, atoms, values)]))
        
        

        class_width = 27
        
        output.append(f"""{"CLASS":<{class_width}} COUNT    MEAN    STDEV""")
        
        divider_idx = len(output)
        output.append("")
        
        first = True
        for display_name, subgroup_data in groups_to_display:
            if not first:
                output.append("")
            first = False
            
            total_count = sum(len(atoms) for _, atoms, _ in subgroup_data)
            
            if len(subgroup_data) > 1:
                output.append(f"{display_name:<{class_width}} {total_count:>5}")
                for cls, atoms, values in subgroup_data:
                    count = len(atoms)
                    mean = np.mean(values) if values else 0
                    stdev = np.std(values) if len(values) > 1 else 0
                    
                    cls_display = self.format_subclass_name(cls, display_name)
                    full_name = f"  - {cls_display}"
                    
                    if values:
                        output.append(f"{full_name:<{class_width}} {count:>5} {mean:>8.3f} {stdev:>7.3f}")
                    else:
                        output.append(f"{full_name:<{class_width}} {count:>5}")
            else:
                cls, atoms, values = subgroup_data[0]
                count = len(atoms)
                mean = np.mean(values) if values else 0
                stdev = np.std(values) if len(values) > 1 else 0
                
                if values:
                    output.append(f"{display_name:<{class_width}} {count:>5} {mean:>8.3f} {stdev:>7.3f}")
                else:
                    output.append(f"{display_name:<{class_width}} {count:>5}")
       
        max_len = 0
        for i, line in enumerate(output):
            if i != divider_idx and not line.startswith("=") and line.strip():
                max_len = max(max_len, len(line))
         
        divider_len = max(max_len + 1, 50)
        output[divider_idx] = "-" * divider_len 
        
        return "\n".join(output)
    
    
    
    def generate_extreme_bonds_table(self):
        """ Generate a table ranking the highest and lowest bond orders for the molecule. """
        
        atom_class = self.parameters.get("atom_class")
        if not atom_class:
            return "No atom classification data available"
        
        bond_params = [key for key in self.parameters.keys() if key.endswith("_bond_orders")]
        
        if not bond_params:
            return "No bond data available"
        
        all_bonds = []
        
        for param_name in bond_params:
            prefix = param_name.replace("_bond_orders", "")
            base_type = self.get_bond_type_display(prefix)
            if base_type is None:
                continue
            
            bond_orders = self.parameters.get(param_name, {})
            bond_lengths_key = param_name.replace("_bond_orders", "_bond_lengths")
            bond_lengths = self.parameters.get(bond_lengths_key, {})
            
            if prefix in ["mbridtm_O", "hetbridm_O", "hetbridtm_O"]:                # Special handling for M-TM Bridging Oxygen (mbridtm_O), which we need to split into (M-O)-TM and M-(O-TM).
                for bond in bond_lengths.keys():
                    order = bond_orders.get(bond, 0.0)
                    length = bond_lengths.get(bond, 0.0)
                    
                    if "-" in bond:
                        atom1, atom2 = bond.split("-")
                        
                        cls1 = atom_class.get(atom1, "")
                        cls2 = atom_class.get(atom2, "")
                        
                        if "Oxygen" in cls1:                                        # Determine which atom oxygen is.
                            other_cls = cls2
                        elif "Oxygen" in cls2:
                            other_cls = cls1
                        else:
                            continue
                        
                        if prefix == "mbridtm_O":
                            if "Framework Metal" in other_cls:                      # Other atom is M, so this is (M-O)-TM.          
                                bond_type = "TM-(O-M)"
                            elif "TM-Substitute" in other_cls:                      # Other atom is TM, so this is M-(O-TM).
                                bond_type = "M-(O-TM)"
                            else:
                                continue
                        elif prefix == "hetbridm_O":
                            if "Framework Metal" in other_cls:
                                bond_type = "X-(O-M)"
                            elif "Heteroatom" in other_cls:
                                bond_type = "M-(O-X)"
                            else:
                                continue
                        elif prefix == "hetbridtm_O":
                            if "Heteroatom" in other_cls:
                                bond_type = "TM-(O-X)"
                            elif "TM-Substitute" in other_cls:
                                bond_type = "X-(O-TM)"
                            else:
                                continue
                    
                        all_bonds.append((bond, order, length, bond_type))
                        
            else:
                for bond in bond_lengths.keys():
                    order = bond_orders.get(bond, 0.0)
                    length = bond_lengths.get(bond, 0.0)
                    all_bonds.append((bond, order, length, base_type))
            
        if not all_bonds:
            return "No bond data available"
            
        sorted_by_order = sorted(all_bonds, key=lambda x: x[1], reverse=True)
            
        highest = sorted_by_order[:10]                                              # Get 10 biggest bond orders
        lowest = sorted_by_order[-10:][::-1]                                        # Get 10 smallest bond orders
            
        
        output = []
        output.append("=" * 63)
        output.append("HIGHEST ORDER                     LOWEST ORDER")
        output.append("=" * 63)
        output.append("")
            
        output.append(f"""{"BOND":<10} {"ORDER":<8} {"TYPE":<10}    {"BOND":<10} {"ORDER":<8} {"TYPE":<20}""")
        output.append("-" * 63)
            
        max_rows = max(len(highest), len(lowest))
            
        for i in range(max_rows):
            if i < len(highest):
                bond, order, length, bond_type = highest[i]
                left = f"{bond:<9} {order:>6.3f}    {bond_type:<9}"
            else:
                left = " " * 32
                
            if i < len(lowest):
                bond, order, length, bond_type = lowest[i]
                right = f"{bond:<9} {order:>6.3f}    {bond_type:<9}"
            else:
                right = " " * 32
                
            output.append(f"{left}     {right}")
            
        return "\n".join(output)
    


    def generate_extreme_charges_table(self, charge_type="qa"):
        """ Generate a table showing the most positive (cations) and most negative (anions) charges.
        Returns formatted table with top 10 most positive and most negative charges. """
        
        if charge_type == "qa":
            charge_values = self.parameters.get("qa_charges")
            charge_label = "QA"
            gap = 16
        #elif charge_type == "hir":
        #    charge_values = self.parameters.get("hir_charges")
        #    charge_label = "Hirshfeld"
        #    gap = 14
        else:
            return f"Unknown charge type: {charge_type}"
        
        atom_class = self.parameters.get("atom_class")
        
        if not atom_class or not charge_values:
            return f"No {charge_label} charge data available"
        
        # Build list of (atom, charge, classification)
        atom_data = []
        for atom, charge in charge_values.items():
            if atom in atom_class:
                classification = atom_class[atom]
                if charge is not None:
                    atom_data.append((atom, charge, classification))
        
        if not atom_data:
            return f"No valid {charge_label} charge data available"
        
        # Sort by charge (most positive to most negative)
        sorted_by_charge = sorted(atom_data, key=lambda x: x[1], reverse=True)
        
        # Get top 10 positive (cations) and top 10 negative (anions)
        cations = sorted_by_charge[:10]
        anions = sorted_by_charge[-10:][::-1]  # Reverse to show most negative first
        
        # Build output
        output = []
        
        # Header with proper spacing - adjust width as needed
        total_width = 100
    
        # Center the header
        header_text_left = f"TOP 10 CATIONS (MOST POSITIVE {charge_label})"
        header_text_right = f"TOP 10 ANIONS (MOST NEGATIVE {charge_label})"
        
        # Build the header line without leading spaces
        header_line = f"{header_text_left}{' ' * gap}{header_text_right}"
        
        # Top border
        output.append("=" * total_width)
        output.append(header_line)
        output.append("=" * total_width)
        output.append("")
        
        # Column headers
        rank_width = 5
        atom_width = 6
        charge_width = 6
        class_width = 25
        
        # Left column: RANK, ATOM, CHARGE, CLASS
        # Right column: RANK, ATOM, CHARGE, CLASS
        left_header = f"{'RANK':<{rank_width}} {'ATOM':<{atom_width}} {'CHARGE':>{charge_width}}    {'CLASS':<{class_width}}"
        right_header = f"{'RANK':<{rank_width}} {'ATOM':<{atom_width}} {'CHARGE':>{charge_width}}    {'CLASS':<{class_width}}"
        
        # Calculate total width of each half
        half_width = rank_width + atom_width + charge_width + class_width + 4
        output.append(f"{left_header:<{half_width}} {right_header:<{half_width}}")
        output.append("-" * total_width)
        
        # Fill rows
        max_rows = max(len(cations), len(anions))
        for i in range(max_rows):
            if i < len(cations):
                atom, charge, cls = cations[i]
                # Format the class name
                cls_display = cls[:class_width-3] + "..." if len(cls) > class_width else cls
                left_row = f"{i+1:<{rank_width}} {atom:<{atom_width}} {charge:>{charge_width}.3f}    {cls_display:<{class_width}}"
            else:
                left_row = " " * (rank_width + atom_width + charge_width + class_width + 3)
            
            if i < len(anions):
                atom, charge, cls = anions[i]
                # Format the class name
                cls_display = cls[:class_width-3] + "..." if len(cls) > class_width else cls
                right_row = f"{i+1:<{rank_width}} {atom:<{atom_width}} {charge:>{charge_width}.3f}    {cls_display:<{class_width}}"
            else:
                right_row = " " * (rank_width + atom_width + charge_width + class_width + 3)
            
            output.append(f"{left_row:<{half_width}} {right_row}")
        
        output.append("=" * total_width)
        
        return "\n".join(output)



    def generate_extreme_spins_table(self, spin_type="loewdin"):
        """ Generate a table showing the largest positive spins and corresponding orbital details. """
        
        if spin_type == "loewdin":
            spin_values = self.parameters.get("loew_spin_pop")
            spin_label = "Löwdin Spin"
        else:
            return f"Unknown spin type: {spin_type}"
        
        atom_class = self.parameters.get("atom_class")
        orbital_data = self.parameters.get("loew_orb_spin")
        
        if not atom_class or not spin_values:
            return f"No {spin_label} spin data available"
            
        atom_data = []
        for atom, spin in spin_values.items():
            if atom in atom_class:
                classification = atom_class[atom]
                if spin is not None:
                    atom_data.append((atom, spin, classification))
        
        if not atom_data:
            return f"No valid {spin_label} spin data available"
            
        sorted_by_spin = sorted(atom_data, key=lambda x: x[1], reverse=True)
        positive = sorted_by_spin[:10]
            
        top_atom_orbital_lines = []
        if orbital_data and positive:
            top_atom = positive[0][0]
            if top_atom in orbital_data:
                data = orbital_data[top_atom]
                try:
                    total_spin = data["total_orbital_spin"]
                    orbitals = data["orbitals"]
                        
                    atom_name = top_atom
                    first_line = f"{atom_name}: Total Spin = {total_spin:.5f}"
                    top_atom_orbital_lines.append(first_line)
                    top_atom_orbital_lines.append("")
                    
                    colon_position = len(atom_name) + 2  
                    sorted_orbitals = sorted(orbitals.items(), key=lambda x: abs(x[1]), reverse=True)

                    for orbital, spin in sorted_orbitals:
                        orbital_str = orbital
                        spaces_needed = colon_position - len(orbital_str)
                        if spaces_needed > 0:
                            orbital_padded = " " * spaces_needed + orbital_str
                        else:
                            orbital_padded = orbital_str
                            
                        top_atom_orbital_lines.append(f"{orbital_padded} : {spin:>10.5f}")
                except (KeyError, TypeError):
                    pass
    
        # Original output building (unchanged)
        output = []
        output.append("=" * 85)
        output.append(f"""{"TOP 10 MOST POSITIVE SPINS":<35} {"TOP ATOM LÖWDIN ORBITAL SPIN DETAILS":>45}""")
        output.append("=" * 85)
        output.append("")
            
        output.append(f"""{"RANK":<6} {"ATOM":<6} {"SPIN":<6} {"CLASS":<20}    {"ORBITAL DETAILS"}""")
        output.append("-" * 85)
            
        max_rows = max(len(positive), len(top_atom_orbital_lines))
        for i in range(max_rows):
            if i < len(positive):
                atom, spin, cls = positive[i]
                cls_display = cls
                if len(cls_display) > 20:
                    cls_display = cls_display[:17] + "..."
                    
                spin_str = f"{spin:>5.3f}"
                left = (f"{i+1:<6} {atom:<6} {spin_str}  {cls_display:<20}")
            else:
                left = " " * 39
        
            if i < len(top_atom_orbital_lines):
                right = top_atom_orbital_lines[i]
            else:
                right = ""
                
            output.append(f"{left:<39}    {right}")
    
        return "\n".join(output)
   


    def generate_spin_summary(self, spin_type="loewdin"):
        """ Generate a spin summary table organized by classification groups. Shows count, mean, and standard deviation for spin populations. """
        
        if spin_type == "loewdin":
            spin_values = self.parameters.get("loew_spin_pop")
            spin_label = "Löwdin Spin Populations"
        else:
            return f"Unknown spin type: {spin_type}"
        
        atom_class = self.parameters.get("atom_class")
        
        if not atom_class or not spin_values:
            return "No spin data available"
        
        output = []
        output.append("=" * 40)
        output.append(f"SPIN SUMMARY ({spin_label})")
        output.append("=" * 40)
        output.append("")
        
        class_to_atoms = {}
        for atom, classification in atom_class.items():
            if classification not in class_to_atoms:
                class_to_atoms[classification] = []
            class_to_atoms[classification].append(atom)
        
        groups_to_display = []
        displayed_classes = set()
        
        for group_name, group_info in class_groups.items():
            classifications = group_info["classifications"]
            subgroup_data = []
            
            for cls in classifications:
                if cls in class_to_atoms:
                    atoms = class_to_atoms[cls]
                    # Get spin values for these atoms
                    values = []
                    for atom in atoms:
                        if atom in spin_values:
                            values.append(spin_values[atom])
                    subgroup_data.append((cls, atoms, values))
                    displayed_classes.add(cls)
            
            if not subgroup_data:
                continue
            
            if "OVERALL" in group_name:
                include = True
            else:
                is_subset_of_overall = False
                for other_name, other_info in class_groups.items():
                    if "OVERALL" in other_name and other_name != group_name:
                        other_classifications = set(other_info["classifications"])
                        if all(cls in other_classifications for cls in classifications):
                            is_subset_of_overall = True
                            break
                include = not is_subset_of_overall
            
            if include:
                display_name = group_name.replace(" - OVERALL", "")
                groups_to_display.append((display_name, subgroup_data))
        
        for classification, atoms in class_to_atoms.items():
            if classification not in displayed_classes:
                values = [spin_values[atom] for atom in atoms if atom in spin_values]
                groups_to_display.append((classification, [(classification, atoms, values)]))
        

        class_width = 27
        
        output.append(f"""{"CLASS":<{class_width}} COUNT    MEAN    STDEV""")
        
        divider_idx = len(output)
        output.append("")
        
        first = True
        for display_name, subgroup_data in groups_to_display:
            if not first:
                output.append("")
            first = False
            
            total_count = sum(len(atoms) for _, atoms, _ in subgroup_data)
            
            if len(subgroup_data) > 1:
                output.append(f"{display_name:<{class_width}} {total_count:>5}")
                for cls, atoms, values in subgroup_data:
                    count = len(atoms)
                    mean = np.mean(values) if values else 0
                    stdev = np.std(values) if len(values) > 1 else 0
                    
                    cls_display = self.format_subclass_name(cls, display_name)
                    full_name = f"  - {cls_display}"
                    
                    if values:
                        output.append(f"{full_name:<{class_width}} {count:>5} {mean:>8.3f} {stdev:>7.3f}")
                    else:
                        output.append(f"{full_name:<{class_width}} {count:>5}")
            else:
                cls, atoms, values = subgroup_data[0]
                count = len(atoms)
                mean = np.mean(values) if values else 0
                stdev = np.std(values) if len(values) > 1 else 0
                
                if values:
                    output.append(f"{display_name:<{class_width}} {count:>5} {mean:>8.3f} {stdev:>7.3f}")
                else:
                    output.append(f"{display_name:<{class_width}} {count:>5}")
        
        max_len = 0
        for i, line in enumerate(output):
            if i != divider_idx and not line.startswith("=") and line.strip():
                max_len = max(max_len, len(line))
        
        divider_len = max(max_len + 1, 50)
        output[divider_idx] = "-" * divider_len
        
        return "\n".join(output)



    def get_bond_type_display(self, prefix):
        """ Convert a parameter prefix to a user-friendly bond type display name. """
        
        # Map prefixes to display names - only oxygen-based ones
        prefix_map = {"hetbridm_O": "X-O-M", "hetbridtm_O": "TM-O-X", "mbridm_O": "M-(O-M)", "mbridtm_O": "M-O-TM", "term_O": "M=O", "termTM_O": "TM=O", "tmtolig_O": "Lig-O-TM",
            # Skip overall groups (return None to skip)
            "all_O": None, "all_M": None, "reg_M": None, "mbridtm_M": None, "het": None, "lig": None, "tmsub": None,}
        
        return prefix_map.get(prefix, prefix)



    def get_parameter(self, param_name):
        """ Retrieve a stored parameter by name. """
        
        if param_name == "all":                                 # Calls all available parameters instead of requested.
            return self.parameters.copy()
        elif param_name == "all_requested":                     # Calls all previously requested parameters.
            return {k: v for k, v in self.parameters.items() if k in self.requested_parameters}
        
        value = self.parameters.get(param_name)
        if value is not None:
            self.requested_parameters.add(param_name)           # Mark called parameter as "requested".
        return value    
    


    def save_to_raw_file(self, filename):
        """ Save all parameters to a text file, writing a comprehensive report that includes:
            - Basic structure data (coordinates, bonds)
            - Electronic properties (QA, Hirshfeld, frontier orbitals)
            - Energy data (electronic, thermal, vibrational)
            - Classification properties for each group """
            
        comp_filename = f"{filename}-Complete.txt"
        
        try:
            with open(comp_filename, "w") as f:

                f.write("\n" + "=" * 50)
                f.write("\nBASIC STRUCTURE DATA\n")
                f.write("=" * 50 + "\n")
                
                if "atom_class" in self.parameters:
                    f.write(self.format_atom_class(self.parameters["atom_class"]))
                    f.write("\n\n\n\n")
                
                if "coords" in self.parameters:
                    f.write("COORDINATES:\n")
                    f.write("-" * 40 + "\n")
                    f.write(f"""{self.format_parameter("coords", self.parameters["coords"])}\n\n""")
                    f.write("\n\n")
                        
                if "bond_lengths" in self.parameters and "bond_orders" in self.parameters:
                    f.write("BOND LENGTHS AND ORDERS:\n")
                    f.write("-" * 45 + "\n")
                    bond_table = self.format_bond_table(self.parameters["bond_lengths"], self.parameters["bond_orders"])
                    f.write(f"{bond_table}\n\n")
                    f.write("\n\n") 
                 
                    
                f.write("\n" + "=" * 50)
                f.write("\nELECTRONIC PROPERTIES\n")
                f.write("=" * 50 + "\n\n")
                
                if "qa_charges" in self.parameters:
                    f.write("\nQA CHARGES:\n")
                    f.write("-" * 50 + "\n")
                    if "sort_qa" in self.parameters:
                        qa_table = self.format_qa_table(self.parameters["qa_charges"], self.parameters["sort_qa"])
                    else:
                        qa_table = self.format_parameter("qa_charges", self.parameters["qa_charges"])
                    f.write(f"{qa_table}\n\n")
                    f.write("\n\n")
                    
                if "hir_charges" in self.parameters and "sort_hir" in self.parameters:
                    f.write("HIRSHFELD CHARGES:\n")
                    f.write("-" * 50 + "\n")
                    if "sort_hir" in self.parameters:
                        hir_table = self.format_hir_table(self.parameters["hir_charges"], self.parameters["sort_hir"])
                    else:
                        hir_table = self.format_parameter("hir_charges", self.parameters["hir_charges"])
                    f.write(f"{hir_table}\n\n")
                    f.write("\n\n")  
                    
                if "loew_spin_pop" in self.parameters and "sort_loew_spin_pop" in self.parameters:
                    f.write("LÖWDIN SPIN POPULATIONS:\n")
                    f.write("-" * 40 + "\n")
                    loewdin_table = self.format_loew_table(self.parameters["loew_spin_pop"], self.parameters["sort_loew_spin_pop"])
                    f.write(f"{loewdin_table}\n\n")
                    f.write("\n\n")    
                    
                if "frontier" in self.parameters and "extended_frontier" in self.parameters:
                    frontier_table = self.format_frontier_data(self.parameters["frontier"], self.parameters["extended_frontier"])
                    f.write(f"{frontier_table}\n\n")
                    f.write("\n\n") 
                    
                if "energy_data" in self.parameters:
                    f.write("ENERGY DATA:\n")
                    f.write("-" * 40 + "\n\n")  
                    energy_table = self.format_elec_energy_data(self.parameters["energy_data"])
                    f.write(f"{energy_table}\n\n")
                    f.write("\n\n")
                    
                if "vib_freq" in self.parameters:
                    f.write("VIBRATIONAL FREQUENCIES:\n")
                    f.write("-" * 40 + "\n\n")
                    vib_display = self.format_vib_freq(self.parameters["vib_freq"])
                    f.write(f"{vib_display}\n\n")
                    f.write("\n\n")
                    

                f.write("\n" + "=" * 50)
                f.write("\nCLASSIFICATION PROPERTIES\n")
                f.write("=" * 50 + "\n\n")
        
                for group_name in class_order:                                  # Iterate through the order defined in Constants.classification_order.
                    prefix = class_to_param.get(group_name)
                    
                    if not prefix:
                        print(f"""Warning: No parameter prefix found for "{group_name}" """)
                        continue
                
                    qa_charges_param = f"{prefix}_qa_charges"
                    qa_sort_param = f"{prefix}_qa_sort"
                    hir_values_param = f"{prefix}_hir_charges"
                    hir_sort_param = f"{prefix}_hir_sort"
                    loew_values_param = f"{prefix}_loew_values"
                    loew_sort_param = f"{prefix}_loew_sort"
                    bond_lengths_param = f"{prefix}_bond_lengths"
                    bond_orders_param = f"{prefix}_bond_orders"
                    
                    has_qa = qa_charges_param in self.parameters
                    has_hir = hir_values_param in self.parameters
                    has_loew = loew_values_param in self.parameters
                    has_bonds = bond_lengths_param in self.parameters
                    
                    if not has_qa and not has_hir and not has_bonds:
                        continue
                    
                    f.write(f"\n{group_name}\n")
                    f.write("=" * 50 + "\n\n")
                    
                    if has_qa:
                        if qa_sort_param in self.parameters:
                            qa_table = self.format_qa_table(self.parameters[qa_charges_param], self.parameters[qa_sort_param])
                        else:
                            qa_table = self.format_qa_charges(self.parameters[qa_charges_param], is_sorted=False)
                        f.write(f"{qa_table}\n\n")
                        
                        values = list(self.parameters[qa_charges_param].values())
                        if values and all(isinstance(v, (int, float)) for v in values):
                            f.write(f"QA Mean: {np.mean(values):.3f}\n")
                            f.write(f"QA Stdev: {np.std(values):.3f}\n")
                        f.write("\n\n")
                    
                    if has_hir:
                        if hir_sort_param in self.parameters:
                            hir_table = self.format_hir_table(self.parameters[hir_values_param], self.parameters[hir_sort_param])
                        else:
                            hir_table = self.format_qa_charges(self.parameters[hir_values_param], is_sorted=False)
                        f.write(f"{hir_table}\n\n")
                        
                        values = list(self.parameters[hir_values_param].values())
                        if values and all(isinstance(v, (int, float)) for v in values):
                            f.write(f"Hir Mean: {np.mean(values):.3f}\n")
                            f.write(f"Hir Stdev: {np.std(values):.3f}\n")
                        f.write("\n\n")
                    
                    if has_loew:
                        if loew_sort_param in self.parameters:
                            loew_table = self.format_loew_table(self.parameters[loew_values_param], self.parameters[loew_sort_param])
                        else:
                            loew_table = self.format_qa_charges(self.parameters[loew_values_param], is_sorted=False)
                        f.write(f"{loew_table}\n\n")
                            
                        values = list(self.parameters[loew_values_param].values())
                        if values and all(isinstance(v, (int, float)) for v in values):
                            f.write(f"Loew Mean: {np.mean(values):.3f}\n")
                            f.write(f"Loew Stdev: {np.std(values):.3f}\n")
                        f.write("\n")
                        
                    if has_bonds and bond_orders_param in self.parameters:
                        bond_table = self.format_bond_table(self.parameters[bond_lengths_param], self.parameters[bond_orders_param])
                        if bond_table and "No bond data" not in bond_table:
                            f.write(f"{bond_table}\n\n")
                    
                    f.write("\n")
                            
            print(f"Parameters saved to {comp_filename}")
            return True
                            
        except Exception as e:
            print(f"Error saving to file: {e}")
            import traceback
            traceback.print_exc()
            return False
     
 

    def save_to_trimmed_file(self, filename):
        """ Save select parameters to a text file, writing a concise report that includes:
            - Basic summary (formula, energy, frontier gap)
            - Atom classifications
            - Charge summaries (QA and Hirshfeld)
            - Extreme charges (most positive/negative)
            - Spin summaries and extreme spins
            - Bond summaries and extreme bonds """
        
        if "atom_class" not in self.parameters:
            print("No atom classification data available. Trimmed file not created")
            return False
        
        trim_filename = f"{filename}-Trimmed.txt"
        
        try:
            with open(trim_filename, "w") as f:
                f.write(self.generate_basic_summary())
                f.write(self.format_atom_class(self.parameters["atom_class"]))
                f.write("\n\n\n\n\n")
                
                #if "hir_charges" in self.parameters and self.parameters["hir_charges"]:
                #    f.write(self.generate_charge_summary("hir"))
                #    f.write("\n\n\n")
                #    f.write(self.generate_extreme_charges_table("hir"))
                #    f.write("\n\n\n\n\n")
                    
                if "qa_charges" in self.parameters and self.parameters["qa_charges"]:
                    f.write(self.generate_charge_summary("qa"))
                    f.write("\n\n\n")
                    f.write(self.generate_extreme_charges_table("qa"))
                    f.write("\n\n\n\n\n")
                    
                if "loew_spin_pop" in self.parameters and self.parameters["loew_spin_pop"]:
                    f.write(self.generate_spin_summary("loewdin"))
                    f.write("\n\n\n")
                    f.write(self.generate_extreme_spins_table("loewdin"))
                    f.write("\n\n\n\n\n")
                    
                f.write(self.generate_bond_summary())
                f.write("\n\n")
                f.write(self.generate_extreme_bonds_table())
                f.write("\n\n\n\n\n")
                
            print(f"Trimmed data saved to {trim_filename}")
            return True
        except Exception as e:
            print(f"Error saving trimmed file: {e}")
            return False
        
        

    def store_parameter(self, param_name, value):
        """ Store a parameter value. """
        
        self.parameters[param_name] = value                     # assigns parameter values with specific names to allow for easier calling
