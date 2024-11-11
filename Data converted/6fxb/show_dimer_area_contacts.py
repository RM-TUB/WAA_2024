import os
import matplotlib.pyplot as plt
import numpy as np
import csv
from matplotlib.ticker import FuncFormatter
from scipy.interpolate import interp1d

def format_ticks(x, pos):
    return '{:.1f}'.format(x)  # Zeit in ns mit einer Dezimalstelle

def parse_csv(file_path):
    """Parses the .csv file and returns time and contact data."""
    time = []
    contacts = []
    
    with open(file_path, 'r') as f:
        csv_reader = csv.reader(f, delimiter=';')
        
        # Skip metadata comments starting with #
        for row in csv_reader:
            if row[0].startswith('#'):
                continue

            # Parse the actual data (assuming first column is time, second is contact data)
            try:
                time.append(float(row[0]))  # Zeit in ns
                contacts.append(float(row[1]))  # Kontakt-Daten
            except ValueError:
                continue
                
    return time, contacts

def find_contact_files(root_directory):
    """Finds relevant .csv files with 'contacts' and '_dimer_area' in the filename."""
    contact_files = []
    for subdir, _, files in os.walk(root_directory):
        for file in files:
            if ('contacts' in file and '_dimer_area' in file and 
                '_inner_dimer_area' not in file and '_outer_dimer_area' not in file and file.endswith('.csv')):
                full_path = os.path.join(subdir, file)
                contact_files.append(full_path)
    return contact_files

def group_files_by_ligand(contact_files):
    """Groups contact files by ligand name."""
    ligands = ['_ca_', '_cab_', '_cad_', '_cla_', '_clab_', '_clad_']
    grouped_files = {ligand: [] for ligand in ligands}

    for file_path in contact_files:
        for ligand in ligands:
            if ligand in file_path:
                grouped_files[ligand].append(file_path)
                break

    return grouped_files

def create_heatmap_for_ligand(ligand, contact_files, root_directory):
    """Creates a heatmap for a specific ligand based on the contact data from the .csv files."""
    all_times = []
    all_contacts = []
    y_labels = []

    for file_path in contact_files:
        time, contacts = parse_csv(file_path)
        
        # Use the filename (without extension) as the Y-axis label
        file_name = os.path.basename(file_path).replace('.csv', '')
        y_labels.append(file_name)

        all_times.append(time)
        all_contacts.append(contacts)

    # Interpolation für höhere Auflösung
    max_time = max([max(t) for t in all_times])
    interpolated_times = np.arange(0, max_time, 0.1)  # 0.1 ns Auflösung
    heatmap_data = []

    for contacts in all_contacts:
        f = interp1d(all_times[0], contacts, bounds_error=False, fill_value=0)  # Interpolation
        interpolated_contacts = f(interpolated_times)
        heatmap_data.append(interpolated_contacts)

    heatmap_data = np.array(heatmap_data)

    # Create the plot
    plt.figure(figsize=(12, len(y_labels) * 0.2))  # Breitere Grafik für bessere Sichtbarkeit
    plt.imshow(heatmap_data, aspect='auto', cmap='gray_r', vmin=0, vmax=1, interpolation='bilinear')  # Verwende bilineare Interpolation
    
    # Add colorbar
    plt.colorbar(label='Number of Contacts')

    # Set labels and ticks
    y_ticks = np.arange(len(y_labels))
    plt.yticks(y_ticks + 0.25, y_labels, fontsize=8)

    # Set X-Ticks und -Labels
    x_tick_indices = np.arange(0, len(interpolated_times), len(interpolated_times) // 10)  # Mehr X-Ticks
    x_tick_labels = np.round(interpolated_times[x_tick_indices], 1)  # Eine Dezimalstelle

    plt.xticks(x_tick_indices, x_tick_labels, fontsize=8)

    # Show the heatmap
    plt.title(f'Contact Heatmap for {ligand} (0 = White, ≥1 = Black)')
    plt.tight_layout()

    # Save the heatmap as a file
    plt.savefig(os.path.join(root_directory, f'heatmap_{ligand.strip("_")}.png'))
    plt.close()

def create_heatmaps(root_directory):
    """Creates heatmaps for each ligand based on the contact data from the .csv files."""
    contact_files = find_contact_files(root_directory)
    grouped_files = group_files_by_ligand(contact_files)

    for ligand, files in grouped_files.items():
        if files:
            print(f"Creating heatmap for ligand: {ligand.strip('_')}")
            create_heatmap_for_ligand(ligand, files, root_directory)

# Skript ausführen
root_directory = "."  # Verzeichnis, in dem das Skript starten soll (aktuelles Verzeichnis)
create_heatmaps(root_directory)
