import os
import fnmatch
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from io import StringIO

# Arrays für Proteine, Liganden und Analysen
proteins = ["2blg", "6fxb"]
ligands = ["ca", "cab", "cad", "cla", "clab", "clad"]
analyses = ["contacts", "mindist", "mindistres", "number_hbonds"]

# Modul: CSV Parsing
def parse_csv(file_path):
    try:
        # Kommentarzeilen mit '#' überspringen
        # Skip die ersten zwei Zeilen für die Meta-Daten, die keine nützlichen Informationen für pandas sind
        df = pd.read_csv(file_path, sep=';', comment='#', skiprows=4, header=0, quotechar='"')
        
        # Überprüfe die Struktur der geladenen Daten
        print(df.head())
        print(df.columns)
        
        return df
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None


# Modul: Subplots erstellen mit Replikaten
def create_subplots(df_list, x_column, y_columns, plot_type="line", title=None):
    num_plots = len(df_list)
    
    if num_plots == 1:
        fig, axs = plt.subplots(1, 1, figsize=(10, 5))  # Eine einzige Achse
        axs = [axs]  # Mache axs zu einer Liste
    else:
        fig, axs = plt.subplots(1, num_plots, figsize=(15, 5))  # Mehrere Achsen

    for i, df in enumerate(df_list):
        print(df.columns)  # Zeigt die Spaltennamen an

        if plot_type == "line":
            for y_col in y_columns:
                axs[i].plot(df[x_column], df[y_col], label=y_col)
        elif plot_type == "box":
            df[y_columns].plot(kind='box', ax=axs[i])
        elif plot_type == "scatter":
            axs[i].scatter(df[x_column], df[y_columns[0]])
        elif plot_type == "bar":
            df.plot(x=x_column, y=y_columns, kind='bar', ax=axs[i])
        elif plot_type == "stacked_bar":
            df.plot(x=x_column, y=y_columns, kind='bar', stacked=True, ax=axs[i])
        elif plot_type == "grouped_bar":
            df.plot(x=x_column, y=y_columns, kind='bar', ax=axs[i])

        axs[i].set_title(f"Replicate {i+1}")
        axs[i].legend()

    if title:
        fig.suptitle(title)
    
    plt.tight_layout()
    plt.show()
# Funktion zum Filtern der Dateien nach Muster
def find_csv_files(directory, proteins, ligands, analyses):
    found_files = []
    
    # Erzeuge Muster für die Suche
    for protein in proteins:
        for ligand in ligands:
            for analyse in analyses:
                pattern = f"{protein}_{ligand}_*_{analyse}_*.csv"
                
                # Durchsuche das Verzeichnis nach dem Muster
                for root, dirs, files in os.walk(directory):
                    for file in fnmatch.filter(files, pattern):
                        found_files.append(os.path.join(root, file))
    
    return found_files

# Funktion zum Batch-Prozessieren der gefundenen Dateien
def process_csv_files(file_list):
    # Gruppiere Dateien nach Analyse
    grouped_files = {}
    
    for file in file_list:
        # Extrahiere den Analyse-Typ aus dem Dateinamen
        for analyse in analyses:
            if analyse in file:
                if analyse not in grouped_files:
                    grouped_files[analyse] = set()  # Verwende ein Set, um Duplikate zu vermeiden
                grouped_files[analyse].add(file)
    
    # Erstelle separate Plots für jede Analyse
    for analyse, files in grouped_files.items():
        print(f"Creating plots for analysis: {analyse}")
        df_list = []
        
        for file in files:
            df = parse_csv(file)
            if df is not None and not df.empty:
                df_list.append(df)
                
        if df_list:
            x_column = df_list[0].columns[0]  # Erste Spalte als X-Achse
            y_columns = df_list[0].columns[1:]  # Weitere Spalten als Y-Achsen
            
            # Subplots für die Replikate erstellen
            create_subplots(df_list, x_column, y_columns, plot_type="line", title=f"Plot for {analyse}")


# Hauptfunktion
def main():
    # Verzeichnis, in dem die Dateien gesucht werden sollen
    directory = '.'  # Pfad anpassen
    
    # Suche nach CSV-Dateien
    csv_files = find_csv_files(directory, proteins, ligands, analyses)
    
    # Verarbeite die gefundenen Dateien
    if csv_files:
        process_csv_files(csv_files)
    else:
        print("Keine passenden Dateien gefunden.")

if __name__ == "__main__":
    main()
