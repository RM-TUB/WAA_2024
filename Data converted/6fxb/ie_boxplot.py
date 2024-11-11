import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Allgemeine Variablen
path = '.'
data_name = r"^6fxb_.*_ie_prod_mol.*\.csv$"
sub_systems = ['ca_1', 'ca_2', 'ca_3', 'cab_1', 'cab_2', 'cab_3',
               'cad_1', 'cad_2', 'cad_3', 'cla_1', 'cla_2', 'cla_3',
               'clab_1', 'clab_2', 'clab_3', 'clad_1', 'clad_2', 'clad_3']
ligands = ['mol1', 'mol2', 'mol3', 'mol4', 'mol5', 'mol6', 'mol7', 'mol8', 'mol9', 'mol10', 'mol11', 'mol12', 'mol13', 'mol14', 'mol15', 'mol16', 'mol17', 'mol18', 'mol19', 'mol20']
main_systems = ['_ca_', '_cab_', '_cad_', '_cla_', '_clab_', '_clad_']

data = {}

# Funktion zur Überprüfung und Ausgabe der CSV-Daten
def check_and_load_data(file_path, system, sub_system):
    try:
        print(f"Versuche, Datei zu laden: {file_path}")
        # Lese die Datei ein
        df = pd.read_csv(file_path, comment='#', sep=';', skiprows=4)
        print(f"Datei erfolgreich geladen: {file_path}")
        print(f"Erste Zeilen der Datei {file_path}:\n{df.head()}")

        if not df.empty and df.shape[1] >= 2:  # Überprüfen, ob der DataFrame nicht leer ist und mindestens 2 Spalten hat
            if system not in data:
                data[system] = {}
            data[system][sub_system] = df
        else:
            print(f"Die Datei {file_path} enthält keine ausreichenden Daten.")
    except Exception as e:
        print(f"Fehler beim Laden der Datei {file_path}: {e}")

# Durchsuche das Verzeichnis nach Dateien
for root, dirs, files in os.walk(path):
    for file in files:
        if re.match(data_name, file):
            for sub_system in sub_systems:
                if sub_system in file:
                    system = sub_system.split('_')[0]
                    file_path = os.path.join(root, file)
                    check_and_load_data(file_path, system, sub_system)

# Überprüfe die Struktur der geladenen Daten
for system in data:
    print(f"\nÜberprüfung der Daten für {system}:")
    for sub_system in data[system]:
        print(f"  - Unterkategorie {sub_system}: {len(data[system][sub_system])} Zeilen geladen.")

# Boxplot-Generierung
def plot_boxplots(data):
    # Definierte Reihenfolge der Systeme für die Plots und deren Farben
    ordered_systems = {
        'ca': '#0e0e96',
        'cla': '#9e8d34',
        'cab': '#2c9187',
        'clab': '#ab7029',
        'cad': '#206b20',
        'clad': '#9e350e'
    }

    # Berechnung der Anzahl der Subplots pro Reihe
    num_systems = len(ordered_systems)
    num_rows = (num_systems + 1) // 2  # 2 Spalten, also Anzahl der Zeilen entsprechend anpassen

    # Setze die Größe der Figur basierend auf der Anzahl der Hauptsysteme
    fig, axes = plt.subplots(nrows=num_rows, ncols=2, figsize=(10, 5 * num_rows), sharey=True)
    
    # Abstand zwischen den Subplots auf 0.05 setzen
    plt.subplots_adjust(hspace=0.05, wspace=0.05)

    # Flatten the axes array for easier iteration
    axes = axes.flatten()

    for idx, (ax, main_system) in enumerate(zip(axes, ordered_systems.keys())):
        # Sammle die Daten für die Boxplots
        categories = ['Coul-SR:Protein-MOL', 'LJ-SR:Protein-MOL']
        box_data = {cat: [] for cat in categories}

        # Daten für den aktuellen Hauptsystem aus dem geladenen DataFrame extrahieren
        for sub_system in data[main_system]:
            df = data[main_system][sub_system]
            if 'Coul-SR:Protein-MOL' in df.columns:
                box_data['Coul-SR:Protein-MOL'].extend(df['Coul-SR:Protein-MOL'].dropna().tolist())
            if 'LJ-SR:Protein-MOL' in df.columns:
                box_data['LJ-SR:Protein-MOL'].extend(df['LJ-SR:Protein-MOL'].dropna().tolist())

        # Führe die Boxplot-Generierung durch
        bp = ax.boxplot(
            [box_data[cat] for cat in categories],
            showmeans=True,
            boxprops=dict(color=ordered_systems[main_system]),  # Boxfarbe
            medianprops=dict(color='#000000'),  # Medianlinie
            whiskerprops=dict(color=ordered_systems[main_system]),  # Whiskerfarbe
            capprops=dict(color=ordered_systems[main_system]),  # Kapfarbe
            meanprops={'marker': '_', 'markerfacecolor': 'white', 'markeredgecolor': 'black', 'markersize': 12}  # Mittelwert
        )
        
        # Setze die X-Achsen-Beschriftungen
        ax.set_xticks(np.arange(1, len(categories) + 1))
        ax.set_xticklabels(categories)

        # Titel innerhalb des Plots anzeigen
        ax.text(0.5, 0.95, f'{main_system}', transform=ax.transAxes,
                horizontalalignment='center', verticalalignment='top', fontsize=10,
                bbox=dict(facecolor='white', alpha=0.8))

        # X-Achsenstriche und -beschriftungen nur in der letzten Reihe anzeigen
        if idx < len(axes) - 2:  # Nicht in der letzten Reihe
            ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        
        # Y-Achsenstriche und -beschriftungen nur in der ersten Spalte anzeigen
        if idx % 2 != 0:  # Nicht in der ersten Spalte
            ax.tick_params(axis='y', which='both', left=False, labelleft=False)
        else:
            ax.set_ylabel('Amount')  # Y-Achse als "Amount" beschriften

        # Anzeige der Boxplotstatistiken
        for i, key in enumerate(categories):
            if box_data[key]:  # Nur wenn Daten vorhanden sind
                y_max = ax.get_ylim()[1]  # Maximaler Y-Wert des aktuellen Axes
                ax.text(i + 1, 0.8 * y_max,
                        f"N: {len(box_data[key])}",
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=8, color='black', bbox=dict(facecolor='white', alpha=0.5))

        # Gitter entfernen
        ax.grid(False)

    plt.tight_layout()
    plt.show()


import matplotlib.pyplot as plt
import numpy as np

def plot_lineplots_LJ(data, window_size=1000, alpha=0.5, y_limits=(-1200, 50)):
    # Definierte Reihenfolge der Systeme für die Plots und deren Farben
    ordered_systems = {
        'ca': '#0e0e96',
        'cab': '#2c9187',
        'cad': '#206b20',
        'cla': '#9e8d34',
        'clab': '#ab7029',
        'clad': '#9e350e'
    }

    # Farbzuteilung für die Subsysteme
    sub_system_colors = {
        'ca_1': '#0e0e96',
        'ca_2': '#0e7b96',
        'ca_3': '#0ea196',
        'cab_1': '#2c9187',
        'cab_2': '#2cb187',
        'cab_3': '#2cc187',
        'cad_1': '#206b20', 
        'cad_2': '#4ca54b',
        'cad_3': '#6fc67e', 
        'cla_1': '#9e8d34', 
        'cla_2': '#9eac34',
        'cla_3': '#9ec734',
        'clab_1': '#ab7029',
        'clab_2': '#ab8d29',
        'clab_3': '#ab9929',
        'clad_1': '#9e350e',
        'clad_2': '#a34e0f',
        'clad_3': '#b15a0e',
    }

    # Definierte Subplot-Reihenfolge
    subplot_order = ['ca', 'cla', 'cab', 'clab', 'cad', 'clad']
    
    # Berechnung der Anzahl der Subplots pro Reihe
    num_systems = len(subplot_order)
    num_rows = (num_systems + 1) // 2  # 2 Spalten, also Anzahl der Zeilen entsprechend anpassen

    # Setze die Größe der Figur basierend auf der Anzahl der Hauptsysteme
    fig, axes = plt.subplots(nrows=num_rows, ncols=2, figsize=(10, 5 * num_rows))
    
    # Abstand zwischen den Subplots auf 0.05 setzen
    plt.subplots_adjust(hspace=0.05, wspace=0.05)

    # Flatten the axes array for easier iteration
    axes = axes.flatten()

    for idx, main_system in enumerate(subplot_order):
        ax = axes[idx]
        
        # Sammle Daten für die Linienplots
        for sub_system in data[main_system]:
            df = data[main_system][sub_system]
            time = df.iloc[:, 0]  # Angenommen, die Zeit ist in der ersten Spalte
            
            if 'LJ-SR:Protein-MOL' in df.columns:
                # Originale Daten mit Transparenz zeichnen
                ax.plot(time, df['LJ-SR:Protein-MOL'], label=sub_system, color=sub_system_colors.get(sub_system, ordered_systems[main_system]), alpha=alpha)
                
                # Gleitenden Durchschnitt berechnen
                rolling_avg = df['LJ-SR:Protein-MOL'].rolling(window=window_size, min_periods=1).mean()
                # Gleitenden Durchschnitt zeichnen
                ax.plot(time, rolling_avg, label=f'{sub_system} (Running Avg)', color='red', linewidth=2, alpha=alpha)

        # Setze Achsenbeschriftungen
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('Amount')
        ax.set_title(f'{main_system} Lennard Jones Potential')
        ax.legend(title='Sub-System', bbox_to_anchor=(1.05, 1), loc='upper left')

        # Gitter entfernen
        ax.grid(False)

        if y_limits is not None and len(y_limits) == 2:
            ax.set_ylim(y_limits[0], y_limits[1])

    plt.tight_layout()
    plt.show()


# Linienplot-Generierungimport matplotlib.pyplot as plt

def plot_lineplots_Coloum(data, y_limits=(-3500, 50), window_size=1000, alpha=0.5):
    # Definierte Reihenfolge der Systeme für die Plots und deren Farben
    ordered_systems = {
        'ca': '#0e0e96',
        'cab': '#2c9187',
        'cad': '#206b20',
        'cla': '#9e8d34',
        'clab': '#ab7029',
        'clad': '#9e350e'
    }

    # Farbzuteilung für die Subsysteme
    sub_system_colors = {
        'ca_1': '#0e0e96',
        'ca_2': '#0e7b96',
        'ca_3': '#0ea196',
        'cab_1': '#2c9187',
        'cab_2': '#2cb187',
        'cab_3': '#2cc187',
        'cad_1': '#206b20', 
        'cad_2': '#4ca54b',
        'cad_3': '#6fc67e', 
        'cla_1': '#9e8d34', 
        'cla_2': '#9eac34',
        'cla_3': '#9ec734',
        'clab_1': '#ab7029',
        'clab_2': '#ab8d29',
        'clab_3': '#ab9929',
        'clad_1': '#9e350e',
        'clad_2': '#a34e0f',
        'clad_3': '#b15a0e',
    }

    # Definierte Subplot-Reihenfolge
    subplot_order = ['ca', 'cla', 'cab', 'clab', 'cad', 'clad']
    
    # Berechnung der Anzahl der Subplots pro Reihe
    num_systems = len(subplot_order)
    num_rows = (num_systems + 1) // 2  # 2 Spalten, also Anzahl der Zeilen entsprechend anpassen

    # Setze die Größe der Figur basierend auf der Anzahl der Hauptsysteme
    fig, axes = plt.subplots(nrows=num_rows, ncols=2, figsize=(10, 5 * num_rows))
    
    # Abstand zwischen den Subplots auf 0.05 setzen
    plt.subplots_adjust(hspace=0.05, wspace=0.05)

    # Flatten the axes array for easier iteration
    axes = axes.flatten()

    for idx, main_system in enumerate(subplot_order):
        ax = axes[idx]
        
        # Sammle Daten für die Linienplots
        for sub_system in data[main_system]:
            df = data[main_system][sub_system]
            time = df.iloc[:, 0]  # Angenommen, die Zeit ist in der ersten Spalte
            if 'Coul-SR:Protein-MOL' in df.columns:
                # Plot für die Subsystemdaten mit Transparenz
                ax.plot(time, df['Coul-SR:Protein-MOL'], label=sub_system, color=sub_system_colors.get(sub_system, ordered_systems[main_system]), alpha=alpha)

                # Berechnung des gleitenden Durchschnitts
                rolling_avg = df['Coul-SR:Protein-MOL'].rolling(window=window_size).mean()
                
                # Gleitenden Durchschnitt plotten mit Transparenz
                ax.plot(time, rolling_avg, label=f'{sub_system} (Running Avg)', color='red', linewidth=2, alpha=alpha)

        # Setze Achsenbeschriftungen
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('Amount')
        ax.set_title(f'{main_system} Coulomb')
        ax.legend(title='Sub-System', bbox_to_anchor=(1.05, 1), loc='upper left')

        # Manuelle Festlegung der Y-Achsen-Limits, falls angegeben
        if y_limits is not None and len(y_limits) == 2:
            ax.set_ylim(y_limits[0], y_limits[1])

        # Gitter entfernen
        ax.grid(False)

    plt.tight_layout()
    plt.show()

# Beispielaufruf der Funktion mit Transparenz, manuell festgelegten Y-Achsen-Limits und einem gleitenden Durchschnitt
# plot_lineplots_Coloum(data, y_limits=(0, 100), window_size=5, alpha=0.5)  # Hier die gewünschten Werte angeben


# Boxplots generieren
plot_boxplots(data)

# Linienplots erzeugen
plot_lineplots_Coloum(data)
plot_lineplots_LJ(data)
