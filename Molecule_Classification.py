# -*- coding: utf-8 -*-
"""
@author: dmalcolm
"""

from collections import defaultdict
from abc import ABC
from math import sqrt

from Config import class_groups
from Parameter_Storage import OrcaParameterStorage
from Physical_Properties import get_element




# ============ CONFIGURATION ============

class_order = list(class_groups.keys())
class_to_param = {name: info["prefix"] for name, info in class_groups.items()}
group_definitions = {name: info["classifications"] for name, info in class_groups.items()}



# ============ BASE CLASSIFIER ============

class BaseClassifier(ABC):
    """ Abstract base class for molecule classifiers. """
    
    def classify_oxygen(self, oxygen_atom, bonded_atoms, framework_metal, heteroatom_classes):
        """ Classify an oxygen atom based on its bonding pattern. """
        pass
    
    def classify_heteroatom(self, atom, bonding_network, framework_metal):
        """ Classify a heteroatom based on its bonding pattern. """
        pass
    
    def classify_framework_metal(self, atom, bonding_network, heteroatom_classes):
        """ Classify a framework metal atom. """
        pass
    
    def classify_proton(self, atom, bonded_atoms, bonding_network, framework_metal, heteroatom_classes, final_classification):
        """ Classify a proton (H atom) based on what it is bonded to. """
        pass
    
    def get_classification_order(self):
        """ Return the preferred order for displaying classifications. """
        pass
    
    def get_classification_groups(self):
        """ Return mapping of group names to classification types. """
        pass

    def get_framework_metals(self):
        """ Return set of elements that can act as framework metals. """
        pass



# ============ POM CLASSIFIER ============

class POMClassifier(BaseClassifier):
    """ Classifier for Polyoxometalate (POM) molecules. """
    
    def __init__(self):
        self.pom_framework_metals = {"Mo", "W", "V", "Nb", "Ta"}
    
    def get_framework_metals(self):
        return self.pom_framework_metals
    
    def classify_oxygen(self, oxygen_atom, bonded_atoms, framework_metal, heteroatom_classes):
        """ Classify POM oxygen atom based on bonding pattern. """
        
        fw_count = 0
        tm_count = 0
        het_count = 0
        ligand_bond = False
        
        for bonded in bonded_atoms:
            bonded_base = get_element(bonded)
            
            if bonded_base == framework_metal:
                fw_count += 1
            elif heteroatom_classes.get(bonded_base) == "TM-Substitute":
                tm_count += 1
            elif heteroatom_classes.get(bonded_base) == "Heteroatom":
                het_count += 1
            elif bonded_base not in [framework_metal, "O"]:
                ligand_bond = True
        
        # POM-specific classification rules
        if het_count >= 1:
            return "TM-Het Bridging Oxygen" if tm_count >= 1 else "M-Het Bridging Oxygen"
        elif fw_count == 2:
            return "M-M Bridging Oxygen"
        elif fw_count == 1 and tm_count == 1:
            return "M-TM Bridging Oxygen"
        elif fw_count == 1:
            return "Terminal Oxygen"
        elif tm_count == 1 and not ligand_bond:
            return "TM-Terminal Oxygen"
        elif tm_count == 1 and ligand_bond:
            return "TM-Ligand Oxygen"
        elif ligand_bond:
            return "Ligand"
        else:
            return "Unclassified Oxygen"
    
    def classify_heteroatom(self, atom, bonding_network, framework_metal):
        """ Classify POM heteroatom based on bonding to oxygen. """
        
        bonded_atoms = bonding_network.get(atom, [])
        oxygens = []
        fw_bonds = []
        
        for bonded in bonded_atoms:
            if get_element(bonded) == "O":
                oxygens.append(bonded)
                fw_bonds.append(self.count_framework_bonds_to_oxygen(bonded, bonding_network, framework_metal))

        # POM-specific heteroatom rules
        if len(oxygens) >= 3 and sum(1 for c in fw_bonds if c >= 3) >= 3:
            return "Heteroatom"
        elif len(oxygens) >= 4 and sum(1 for c in fw_bonds if c >= 1) >= 4:
            return "TM-Substitute"
        else:
            return "Ligand"
    
    def classify_framework_metal(self, atom, bonding_network, heteroatom_classes):
        """ Classify POM framework metal. """
        # Check for M-O-TM bridging
        for bonded in bonding_network.get(atom, []):
            if get_element(bonded) == "O":
                for oxygen_bond in bonding_network.get(bonded, []):
                    if oxygen_bond != atom:
                        bond_base = get_element(oxygen_bond)
                        if heteroatom_classes.get(bond_base) == "TM-Substitute":
                            return "TM-Bridging Framework Metal"
        return "Regular Framework Metal"
    
    def classify_proton(self, atom, bonded_atoms, bonding_network, framework_metal, heteroatom_classes, final_classification):
        """ Classify POM proton. """
        
        for bonded in bonded_atoms:
            bonded_base = get_element(bonded)
            if bonded_base == "O":
                if bonded in final_classification:
                    oxy_class = final_classification[bonded]
                else:
                    oxy_class = self.classify_oxygen(bonded, bonding_network.get(bonded, []), framework_metal, heteroatom_classes)
                
                if "Ligand" in oxy_class:
                    return "Ligand"
                else:
                    return "Proton"
        return "Ligand"
    
    def count_framework_bonds_to_oxygen(self, oxygen_atom, bonding_network, framework_metal):
        """ Count framework metals bonded to an oxygen. """
        
        count = 0
        for bonded in bonding_network.get(oxygen_atom, []):
            if get_element(bonded) == framework_metal:
                count += 1
        return count
    
    def get_classification_order(self):
        """ Get POM-specific classification display order. """
        
        return class_order
    
    def get_classification_groups(self):
        """ Get POM-specific classification groups. """
        return group_definitions
    

    
# ============ CLASSIFIER REGISTRY ============

classifier_library = {"pom": POMClassifier,}


def get_classifier(molecule_type: str = "pom"):
    """ Get a classifier instance for the specified molecule type. """
    
    if molecule_type not in classifier_library:
        raise ValueError(f"Unknown molecule type: {molecule_type}. "
                        f"Available types: {list(classifier_library.keys())}")
    
    return classifier_library[molecule_type]()


    
# ============ HELPER FUNCTIONS ============
    
def desig_prop(results, parameter_storage=None):
    """ Process designation properties using the order from Constants """
    
    if parameter_storage is None:
        parameter_storage = OrcaParameterStorage()

    for group_name, group_info in class_groups.items():
        designations = group_info["classifications"]
        param_prefix = group_info["prefix"]    
        
        prop_data = filter_properties(results, designations)
        
        bond_lengths, bond_orders, qa_charges, qa_mean, qa_stdev, hir_charges, hir_mean, hir_stdev = prop_data

        parameter_storage.store_parameter(f"{param_prefix}_bond_lengths", bond_lengths)
        parameter_storage.store_parameter(f"{param_prefix}_bond_orders", bond_orders)
        
        parameter_storage.store_parameter(f"{param_prefix}_qa_charges", qa_charges)
        sorted_qa_group = dict(sorted(qa_charges.items(), key=lambda x: x[1]))
        parameter_storage.store_parameter(f"{param_prefix}_qa_sort", sorted_qa_group)
        parameter_storage.store_parameter(f"{param_prefix}_qa_mean", qa_mean)
        parameter_storage.store_parameter(f"{param_prefix}_qa_stdev", qa_stdev)
        
        parameter_storage.store_parameter(f"{param_prefix}_hir_charges", hir_charges)
        sorted_hir_group = dict(sorted(hir_charges.items(), key=lambda x: x[1]))
        parameter_storage.store_parameter(f"{param_prefix}_hir_sort", sorted_hir_group)
        parameter_storage.store_parameter(f"{param_prefix}_hir_mean", hir_mean)
        parameter_storage.store_parameter(f"{param_prefix}_hir_stdev", hir_stdev)

    
    return results



def dict_values_stats(data):
    """ Calculates mean and population standard deviation of all values in a dictionary. """ 
    
    if not data:
        return 0.0, 0.0
    if len(data) == 1:
        return round(next(iter(data.values())), 4), 0.0
    
    values = list(data.values())
    mean = sum(values) / len(values)
    variance = sum((x - mean) ** 2 for x in values) / len(values)
    return round(mean, 4), round(sqrt(variance), 4)



def filter_by_class(data_dict, class_dict, target_class):
    """ Filters dictionary keeping items where any atom matches target_class. 
    Returns dict of filtered items {atom/bond: float, ...}
    e.g. When target_class is "Terminal Oxygen", returns data for atom/bond with corresponding tag. """
    
    result = {}
    missing = set()
  
    for key, value in data_dict.items():
        atoms = key.split("-")
        matched = False
        
        for atom in atoms:
            if atom in class_dict:
                if class_dict[atom] == target_class:
                    matched = True
                    break
            else:
                missing.add(atom)
        
        if matched:
            result[key] = value
    
    return result



def filter_properties(results, atom_class, designations):
    """ Filter properties for given designations. """
            
    data_map = {"bond_l": results["bond_lengths"], "bond_o": results["bond_orders"], "qa": results["qa_charges"], "hir": results["hir_charges"]}
            
    combined = {"bond_l": {}, "bond_o": {}, "qa": {}, "hir": {}}
    atom_class = results["atom_class"]
            
    for designation in designations:
        for key in combined.keys():
            filtered = filter_by_class(data_map[key], atom_class, designation)
            combined[key].update(filtered)
    
    for bond in combined["bond_l"].keys():          # For bonds in combined bond_l that are missing from combined bond_o, add them with default order 0.0.
        if bond not in combined["bond_o"]:
            combined["bond_o"][bond] = 0.0
            
    sorted_data = {k: dict(sorted(v.items())) for k, v in combined.items()}
    qa_mean, qa_stdev = dict_values_stats(combined["qa"])
    hir_mean, hir_stdev = dict_values_stats(combined["hir"])
            
    return sorted_data["bond_l"], sorted_data["bond_o"], sorted_data["qa"], qa_mean, qa_stdev, sorted_data["hir"], hir_mean, hir_stdev
             


def get_atoms_by_class(classifications, target_class):
    """ Get all atoms with a specific classification.
    Returns a list of atoms with the corresponding target_class. """
    
    return [atom for atom, cls in classifications.items() if cls == target_class]

    

# ============ GENERIC CLASSIFICATION FUNCTIONS ============

def build_bonding_network(bond_lengths):
    """ Build bonding network from bond length dictionary. 
    Returns all_atoms dict {"H0", "H1", ...}.
    Returns bonding_network dict that contains atom as a key, with all the other atoms it's bound to as the items {"C5": ["C6", "C7", "C9", "O8"]}. """
    
    all_atoms = set()
    bonding_network = defaultdict(list)
    
    for bond in bond_lengths:
        try:
            a1, a2 = bond.split("-")
            all_atoms.add(a1)
            all_atoms.add(a2)
            bonding_network[a1].append(a2)
            bonding_network[a2].append(a1)
        except ValueError:
            print(f"WARNING: Skipping malformed bond key: {bond}")
            continue
    
    return all_atoms, dict(bonding_network)



def classify_atoms(element_counts, bond_lengths, molecule_type: str = "pom"):
    """ Classify all atoms in a molecule based on current bonding network and classifier. 
    Returns dict {atom: designation, ...}. 
    e.g. {"O0": "Terminal Oxygen", "H1": "Ligand", ...}. """
    
    if not element_counts or not bond_lengths:
        return {}
    
    classifier = get_classifier(molecule_type)
    
    all_atoms, bonding_network = build_bonding_network(bond_lengths)
    
    framework_metals = classifier.get_framework_metals()                            # Determine framework metal (most abundant framework metal).
    candidates = [(e, c) for e, c in element_counts.items() if e in framework_metals]
    
    if not candidates:
        print("WARNING: No framework metals found! Using first metal as default.")
        framework_metal = next(iter(element_counts.keys()))
    else:
        framework_metal = max(candidates, key=lambda x: x[1])[0]
    
    heteroatom_classes = {}                                                         # Classify heteroatoms (non-O, non-framework elements).
    for atom in all_atoms:
        base_element = get_element(atom)
        if base_element not in ("O", framework_metal):
            heteroatom_classes[base_element] = classifier.classify_heteroatom(atom, bonding_network, framework_metal)
    
    final = {}                                                                      # Classify all atoms.
    protonated_oxygens = set()
    
    for atom in all_atoms:
        base_element = get_element(atom)
        bonded = bonding_network.get(atom, [])
        
        if base_element == "H":
            final[atom] = classifier.classify_proton(atom, bonded, bonding_network, framework_metal, heteroatom_classes, final)
        
        elif base_element == "O":
            final[atom] = classifier.classify_oxygen(atom, bonded, framework_metal, heteroatom_classes)
            if any(get_element(b) == "H" for b in bonded):                          # Track protonated oxygens.
                protonated_oxygens.add(atom)
        
        elif base_element == framework_metal:
            final[atom] = classifier.classify_framework_metal(atom, bonding_network, heteroatom_classes)
                    
        else:
            final[atom] = heteroatom_classes.get(base_element, "Ligand")
    
    return dict(sorted(final.items(), key=lambda x: extract_number(x[0])))          # Sort by atom number.



def count_elements(atom_labels):
    """ Count occurrences of each element in atom labels. 
    Return dict {atom: count, ...}. """
    
    counts = {}
    for label in atom_labels:
        e = get_element(label)
        if e:
            counts[e] = counts.get(e, 0) + 1
    return counts



def extract_number(atom_label):
    """ Extract numeric suffix from atom label.
    e.g. Takes 'Fe1', returns integer 1. """
    
    num_str = "".join(filter(str.isdigit, atom_label))
    return int(num_str) if num_str else 0



def format_superscript_charge(charge_val):
    """ Format charge as Unicode superscript. """
        
    superscript_conv = {"0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹"}
        
    abs_charge = abs(charge_val)
    charge_str = "".join(superscript_conv[digit] for digit in str(abs_charge))
    
    if charge_val > 0:
        return charge_str + "⁺"
    else:
        return charge_str + "⁻"



def generate_molecule_formula(atom_class, charge=None):
    """ Generate the molecule formula from the atom classification """

    if not atom_class:
        return None
    
    subscript_conv = {"0": "₀", "1": "₁", "2": "₂", "3": "₃", "4": "₄", "5": "₅", "6": "₆", "7": "₇", "8": "₈", "9": "₉"}

    def format_subscript(number):
        """ Convert a number to Unicode subscript. """
        
        if number <= 1:
            return ""
        return "".join(subscript_conv[digit] for digit in str(number))
    
    hetero_classes = group_definitions.get("HETEROATOM", [])                                                                                                       # Define category -> list of classification strings
    framework_classes = group_definitions.get("FRAMEWORK METAL - OVERALL", [])
    tm_sub_classes = group_definitions.get("TM-SUBSTITUTE", [])
    oxygen_classes = group_definitions.get("OXYGEN - OVERALL", [])
    ligand_classes = (group_definitions.get("LIGAND - OVERALL", []) + group_definitions.get("LIGAND - OXYGEN (Lig-O-TM)", []) + ["Proton"])

    categories = {"hetero": hetero_classes, "framework": framework_classes, "tm_sub": tm_sub_classes, "oxygen": oxygen_classes, "ligand": ligand_classes}          # Counters per category: element -> count

    counts = {}
    assigned = set()

    for category, class_list in categories.items():
        atoms = [atom for atom, cls in atom_class.items() if cls in class_list]
        counts[category] = count_elements(atoms)
        assigned.update(atoms)

    unassigned = [atom for atom in atom_class.keys() if atom not in assigned]                                                                                      # Fallback: any unassigned atoms go to ligand
    if unassigned:
        counts["ligand"] = {**counts.get("ligand", {}), **count_elements(unassigned)}

    def part_string(counter):
        """ Build a string from a counter with Unicode subscripts. """
        
        result = ""
        for elem, count in sorted(counter.items()):
            result += elem
            if count > 1:
                result += format_subscript(count)
        return result

    hetero_part = part_string(counts["hetero"])
    framework_part = part_string(counts["framework"])
    tm_sub_part = part_string(counts["tm_sub"])
    
    total_o = sum(counts["oxygen"].values())
    oxygen_part = f"O{format_subscript(total_o)}" if total_o > 0 else ""
    
    ligand_str = part_string(counts["ligand"])
    ligand_part = f"({ligand_str})" if ligand_str else ""

    formula = hetero_part + framework_part + tm_sub_part + oxygen_part + ligand_part
    
    return formula if formula else None





# =============================================================================
# MAIN EXECUTION (For Testing)
# =============================================================================

if __name__ == "__main__":
    # Test with POM example
    test_atom_labels = ["Mo1", "Mo2", "O1", "O2", "O3", "P1"]
    test_bond_lengths = { "Mo1-O1": 1.72, "Mo1-O2": 1.82, "Mo1-O3": 1.91, "Mo2-O1": 1.73, "Mo2-O2": 1.83, "P1-O3": 1.85}
    
    element_counts = count_elements(test_atom_labels)
    classifications = classify_atoms(element_counts, test_bond_lengths)
    
    print("POM Classification Test:")
    
    summary = {}
    for cls in classifications.values():
        summary[cls] = summary.get(cls, 0) + 1
    summary = dict(sorted(summary.items()))
    
    print("\n" + "=" * 60)
    print("ATOM CLASSIFICATION REPORT")
    print("=" * 60)
    
    # Print summary
    print("\nClassification Summary:")
    for cls, count in summary.items():
        print(f"  {cls}: {count}")
    
    # Print detailed atom list by classification
    print("\nDetailed Atom List:")
    for cls, count in sorted(summary.items()):
        atoms = get_atoms_by_class(classifications, cls)
        atom_list = ", ".join(sorted(atoms, key=extract_number))
        print(f"\n  {cls} ({count}):")
        print(f"    {atom_list}")
    
    print("\n" + "=" * 60)