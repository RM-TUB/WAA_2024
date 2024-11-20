import os
import pandas as pd

parsed_files = [""]

file_names = ["2blg_clab_1", 
    "2blg_clab_2_",
    "2blg_clab_3_",
    "2blg_cla_1_",
    "2blg_cla_2_",
    "2blg_cad_1_",
    "2blg_ca_2_",
    "2blg_ca_2_prod_mindistres_mol1_55_200_21200.csv",
    "2blg_clad_1_",
    "2blg_clad_2_",
    "6fxb_clad_1_",
    "6fxb_clad_2_",
    "6fxb_clad_3_",
    "6fxb_clab_1_",
    "6fxb_clab_2_",
    "6fxb_clab_3_",
    "6fxb_cla_1_",
    "6fxb_cla_3_",
    "6fxb_cad_1_",
    "6fxb_cad_2_",
    "6fxb_cad_3_",
    "6fxb_cab_1_",
    "6fxb_cab_2_",
    "6fxb_cab_3_"]

ligand_numbers = ["mol2", 
    "mol9_",
    "mol1_",
    "mol6_",
    "mol6_",
    "mol8_",
    "mol1_",
    "mol1_",
    "mol8_",
    "mol9_",
    "mol15_",
    "mol15_",
    "mol15_",
    "mol8_",
    "mol20_",
    "mol13_",
    "mol15_",
    "mol2_",
    "mol17_",
    "mol5_",
    "mol3_",
    "mol16_",
    "mol16_",
    "mol16_"]

ligand_orders = ["ca", "cab", "cad", "cla", "clab", "clad"]

protein_orders = ["6fxb", "2blg"]

replicates = ["1", "2", "3"]

output_data =  {file_name: [] for file_name in file_names}

def parse_csv(file_path):
    """Parses the .csv file and returns time and contact data."""
    residues = []
    with open(file_path, 'r') as f:               
        try:

            df = pd.read_csv(file_path, sep=';', skiprows=5, header=None)
            first_column = df.iloc[:, 0].tolist()
            second_column = df.iloc[:, 1].tolist()

            for row0, row1 in zip(first_column, second_column):
                if row1 <= 0.55:
                    residues.append(row0)        
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
    for file_name in file_names:
        if file_name in os.path.basename(file_path):
            output_data[file_name].append(residues)

def find_contact_files():
    """Finds relevant .csv files with 'contacts' and '_dimer_area' in the filename."""
    print("function called")
    root_directory = '.'
    contact_files = []

    # Walk through directories to find matching files
    for subdir, _, files in os.walk(root_directory):
        for file in files:
            for file_name, ligand_number in zip(file_names, ligand_numbers):
                if (file_name in file and ligand_number in file and 'mindistres' in file and file.endswith('.csv') and '_r_' not in file and 'neg' not in file and 'pos' not in file):
                    full_path = os.path.join(subdir, file)
                    contact_files.append(full_path)

    # Sort by protein, ligand, replicate, and mol number
    contact_files.sort(key=lambda x: (
        next((i for i, fn in enumerate(file_names) if fn in x), len(file_names)),  # sort by file name
        next((i for i, ln in enumerate(ligand_numbers) if ln in x), len(ligand_numbers)),  # sort by ligand number
        int(next((rep for rep in replicates if f'_{rep}_' in x), '0'))  # sort by replicate as an integer
    ))
   
    for contact_file in contact_files:             
        if os.path.basename(contact_file) not in parsed_files:
            parse_csv(contact_file)
            parsed_files.append(os.path.basename(contact_file))
            print(f"parsed {contact_file}")

def generate_pymol_select_commands(output_data):
    select_commands = []
    
    for file_name, residues in output_data.items():
        if residues:  # Überprüfen, ob es Residuen gibt
            # Erstelle den Befehl für die Auswahl
            selection = f"select {os.path.basename(file_name)}, (res "
            residue_indices = set()  # Verwenden von set, um Duplikate zu vermeiden
            
            # Alle Residuen in ein Set hinzufügen und in ganze Zahlen umwandeln
            for residue_list in residues:
                residue_indices.update(int(residue) for residue in residue_list)  # alle Residuen zu einer Liste hinzufügen
            
            # Sortieren und in Intervalle umwandeln
            sorted_indices = sorted(residue_indices)
            ranges = []
            start = sorted_indices[0]
            prev = sorted_indices[0]
            
            for idx in sorted_indices[1:]:
                if idx != prev + 1:  # Wenn der Index nicht aufeinander folgt
                    ranges.append(f"{start}-{prev}" if start != prev else str(start))
                    start = idx
                prev = idx
            
            # Letzte Range hinzufügen
            ranges.append(f"{start}-{prev}" if start != prev else str(start))
            
            # Füge die Bereiche zum Auswahlbefehl hinzu
            selection += ','.join(ranges) + ")"
            select_commands.append(selection)
    
    return '\n'.join(select_commands)



find_contact_files()
pymol_select_commands = generate_pymol_select_commands(output_data)
print(pymol_select_commands)