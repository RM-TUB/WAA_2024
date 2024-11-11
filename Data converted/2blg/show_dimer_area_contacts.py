import os
import matplotlib.pyplot as plt
import numpy as np
import csv
from matplotlib.ticker import FuncFormatter

def format_ticks(x, pos):
    return '{:.0f}'.format(x / 1000)  # Umrechnung in Tausend

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
                time.append(float(row[0]) / 1000.0)  # Convert time from ps to ns
                contacts.append(float(row[1]))  # Contact data
            except ValueError:
                # In case of a conversion error, skip the row
                continue
                
    return time, contacts

def find_contact_files(root_directory):
    """Finds relevant .csv files with 'contacts' and '_dimer_area' in the filename, excluding '_inner_dimer_area' and '_outer_dimer_area'."""
    contact_files = []
    for subdir, _, files in os.walk(root_directory):
        for file in files:
            # Check if 'contacts' and '_dimer_area' are in the filename and exclude unwanted terms
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
                break  # File belongs to one ligand, stop checking others

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

    # Determine the maximum number of time points for the x-axis
    max_time_points = max([len(t) for t in all_times])

    # Pad all contact arrays to have the same length
    padded_contacts = []
    for contacts in all_contacts:
        padded_contacts.append(np.pad(contacts, (0, max_time_points - len(contacts)), constant_values=0))

    # Create the heatmap data
    heatmap_data = np.array(padded_contacts)

    # Create the plot
    plt.figure(figsize=(10, len(y_labels) * 0.2))  # Height of the figure adjusted
    plt.imshow(heatmap_data, aspect='auto', cmap='gray_r', vmin=0, vmax=1, interpolation='nearest')  # White at 0, black at ≥1
    
    # Add colorbar
    plt.colorbar(label='Number of Contacts')

    # Set labels and ticks
    y_ticks = np.arange(len(y_labels))  # Y-Ticks entsprechend der Anzahl der Labels
    plt.yticks(y_ticks + 0.25, y_labels, fontsize=8)  # Adjust Y-Ticks for spacing

    # Generiere die X-Werte basierend auf den Zeitdaten
    if all_times:
        x_values = all_times[0]  # Verwende die Zeitdaten des ersten Datensatzes

        # Setze die X-Ticks und -Labels
        num_ticks = 5
        x_tick_indices = np.linspace(0, len(x_values) - 1, num_ticks).astype(int)
        x_tick_labels = np.round(np.array(x_values)[x_tick_indices]).astype(int)  # Zeit in ns ohne Nachkommastellen

        plt.xticks(x_tick_indices, x_tick_labels, fontsize=8)  # X-Achse in ns

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

    # Create a heatmap for each ligand
    for ligand, files in grouped_files.items():
        if files:
            print(f"Creating heatmap for ligand: {ligand.strip('_')}")
            create_heatmap_for_ligand(ligand, files, root_directory)

# Skript ausführen
root_directory = "."  # Verzeichnis, in dem das Skript starten soll (aktuelles Verzeichnis)
create_heatmaps(root_directory)
