import pandas as pd
import os
import re

n_Ligand = 20

# Verzeichnis, in dem die CSV-Dateien gespeichert sind
directory = "."  # Pfad zu deinem Verzeichnis

# Reguläre Ausdrücke für die Dateinamen
data_name_pattern = r"^6fxb_.*_ie_prod_mol.*\.csv$"
sasa_name_pattern = r"^6fxb_.*prod_sasa.*\.csv$"

# Listen für die Dateinamen
data_files = []
sasa_files = []

# Durchlaufe das Verzeichnis und finde passende Dateien
for filename in os.listdir(directory):
    if re.match(data_name_pattern, filename):
        data_files.append(filename)
    elif re.match(sasa_name_pattern, filename):
        sasa_files.append(filename)

# Überprüfe, ob die Dateien gefunden wurden
if not data_files or not sasa_files:
    print("Keine passenden Dateien gefunden.")
else:
    for data_file in data_files:
        # Lade die IE-Daten
        data_df = pd.read_csv(os.path.join(directory, data_file), sep=';', skiprows=4)

        # Überprüfen Sie die Spaltennamen
        print(f"Spaltennamen in {data_file}: {data_df.columns.tolist()}")

        # Extrahiere das Präfix für die SASA-Datei
        prefix = '_'.join(data_file.split('_')[:-3])  # z.B. "6fxb_clad_3"

        # Suche nach dem passenden SASA-Dateinamen
        corresponding_sasa_file = next((sasa for sasa in sasa_files if sasa.startswith(prefix)), None)
        
        if corresponding_sasa_file:
            # Lade die SASA-Daten
            sasa_df = pd.read_csv(os.path.join(directory, corresponding_sasa_file), sep=';', skiprows=4)
            
            # Überprüfen Sie die Spaltennamen der SASA-Daten
            print(f"Spaltennamen in {corresponding_sasa_file}: {sasa_df.columns.tolist()}")

            # Überprüfe, ob die Zeitspalten übereinstimmen
            if 'Time (ps)' in data_df.columns and 'Time (ps)' in sasa_df.columns:
                if not (data_df['Time (ps)'] == sasa_df['Time (ps)']).all():
                    print(f"Zeitspalten stimmen nicht überein für {data_file} und {corresponding_sasa_file}.")
                    continue
            else:
                print(f"Spalte 'Time (ps)' fehlt in einer der Dateien: {data_file} oder {corresponding_sasa_file}.")
                continue

            # Normalisierung: Beide Y-Spalten durch die Total-Spalte aus der SASA-Datei teilen
            normalized_df = data_df.copy()
            normalized_df['Coul-SR:Protein-MOL'] = data_df['Coul-SR:Protein-MOL'] / sasa_df['Total'] / n_Ligand
            normalized_df['LJ-SR:Protein-MOL'] = data_df['LJ-SR:Protein-MOL'] / sasa_df['Total'] / n_Ligand

            # Erstelle den neuen Dateinamen
            new_filename = re.sub(r"(_ie_prod_mol.*\.csv$)", r"_ie_prod_mol_norm\1", data_file)
            # Speichere die normalisierten Daten in einer neuen CSV-Datei
            normalized_df.to_csv(os.path.join(directory, new_filename), sep=';', index=False)
            print(f"Normalisierte Daten gespeichert in: {new_filename}")
        else:
            print(f"Kein passendes SASA-Datei gefunden für {data_file}.")
