import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import colorsys
import re

def extract_replica_number(file_path):
    match = re.search(r'_ca_(\d+)', file_path)  # Beispielregex, um die Zahl nach "_ca_" zu extrahieren
    if match:
        return int(match.group(1))
    return 0


# Funktion zum Parsen von .xvg-Dateien
def parse_xvg(file_path):
    """
    Liest eine .xvg-Datei und gibt Zeit- und Datenwerte zurück.

    :param file_path: Pfad zur .xvg-Datei
    :return: Tuple von Arrays (Zeitwerte, Datenwerte)
    """
    x = []
    y_values = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('@') and not line.startswith('#'):
                data = line.split()
                if len(data) > 1:
                    # Zeit von Pikosekunden auf Nanosekunden konvertieren
                    x.append(float(data[0]))
                    if not y_values:
                        y_values = [[] for _ in range(len(data) - 1)]
                    for i in range(len(data) - 1):
                        y_values[i].append(float(data[i + 1]))
    return np.array(x), np.array(y_values).T

def plot_data(p1_files, p2_files, title, output_file, data_index=0, ylabel='Value [nm]', colors=None, alpha=1, x_label='Time [ns]', legend_labels=None, constrained_layout=True):
    ligands = ['CA_', 'CLA_', 'CAB', 'CLAB', 'CAD', 'CLAD', 'NOLIGAND']
    num_cols = 2
    num_rows = 3
    figsize_x = 8.27
    figsize_y = 11.69
    fig_dpi = 600
    fontsize_ticks = 12
    fontsize_textbox = 12
    fontsize_notes = 8
    fontsize_labels = 12
    fontsize_title = 14

    
    colors = {
        'CA_': '#00FFFF',
        'CAB': '#00CCCC',
        'CAD': '#009999',
        'CLA_': '#FFFF00',
        'CLAB': '#FFCC00',
        'CLAD': '#FF9900',
        'NOLIGAND': '#888888',
    }
    if legend_labels is None:
        legend_labels = {
            'CA_': 'CA',
            'CLA_': 'CLA',
            'CAB': 'CAB',
            'CLAB': 'CLAB',
            'CAD': 'CAD',
            'CLAD': 'CLAD',
            'NOLIGAND' : ''
        }

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(figsize_x, figsize_y), dpi = fig_dpi)

    if colors is None:
        colors = cm.viridis(np.linspace(0, 1, len(ligands)))
    
    global_y_min, global_y_max = float('inf'), float('-inf')
    for file_path in p1_files + p2_files:
        x, y_values = parse_xvg(file_path)
        global_y_min = min(global_y_min, min(y_values[:, data_index]))
        global_y_max = max(global_y_max, max(y_values[:, data_index]))

    if 'RMSF' in title.upper():
        upper_limit = 0.64
        lower_limit = 0
        x_upper_limit = 163

    elif 'RMSD' in title.upper():
        upper_limit = 0.375
        lower_limit = 0
        x_upper_limit = 301
    else:
        y_range = global_y_max - global_y_min
        upper_limit = global_y_min + y_range * 1.5
        upper_limit = 1.565
        lower_limit = 1.465
        x_upper_limit = 301

    y_ticks = np.linspace(lower_limit, upper_limit, 5)

    for i, ligand in enumerate(ligands[:-1]): 
        row = i // num_cols
        col = i % num_cols
        ax = axes[row, col]
        
        p1_files_ligand = [f for f in p1_files if ligand in f]
        p2_files_ligand = [f for f in p2_files if ligand in f]
        
        p1_files_ligand.sort(key=extract_replica_number)
        p2_files_ligand.sort(key=extract_replica_number)

        for idx, p1_file in enumerate(p1_files_ligand):
            alpha = [1.0, 0.5, 0.25]
            x_p1, y_p1 = read_xvg(p1_file)
            label=f'BLG-{legend_labels[ligand]} ({idx+1})'
            graphs = ax.plot(x_p1, y_p1[:, data_index], color=colors[ligand], alpha=alpha[idx], label=label, linewidth=0.5)
            if 'RMSF' in title.upper():
                for graph in graphs:
                    graph.set_linewidth(1)
            
#        for idx, p2_file in enumerate(p2_files_ligand):
#            adjust_alpha = 1 - (((1 + idx) - (1)) / (3 - 0))
#            x_p2, y_p2 = read_xvg(p2_file)
#            label=f'BLG$_{{{2}}}$-{legend_labels[ligand]} ({idx+1})'
#            ax.plot(x_p2, y_p2[:, data_index], color=colors[ligand], alpha=alpha[idx], linestyle='--', label=label)
        

        p1_files_ligand = [f for f in p1_files if 'NOLIGAND' in f]
        p2_files_ligand = [f for f in p2_files if 'NOLIGAND' in f]
        
        p1_files_ligand.sort(key=extract_replica_number)
        p2_files_ligand.sort(key=extract_replica_number)

        for idx, p1_file in enumerate(p1_files_ligand):
            alpha = [0.6, 0.3, 0.15]
            x_p1, y_p1 = read_xvg(p1_file)
            label=f'BLG-NK ({idx+1})'
            graphs = ax.plot(x_p1, y_p1[:, data_index], color=colors['NOLIGAND'], alpha=alpha[idx], label=label, linewidth=0.5)
            if 'RMSF' in title.upper():
                for graph in graphs:
                    graph.set_linewidth(1)        
#        for idx, p2_file in enumerate(p2_files_ligand):
#            adjust_alpha = 1 - (((1 + idx) - (1)) / (3 - 0))
#            x_p2, y_p2 = read_xvg(p2_file)
#            label=f'BLG$_{{{2}}}$-{legend_labels['NOLIGAND']} ({idx+1})'
#            ax.plot(x_p2, y_p2[:, data_index], color=colors['NOLIGAND'], alpha=alpha[idx+3], linestyle='--', label=label)


        if row == num_rows - 1:
            ax.set_xlabel(x_label, fontsize=fontsize_labels)
        else:
            ax.set_xticklabels([])
            ax.tick_params(axis='x', which='both', bottom=False)

        if col == 0 and row == 1:
            ax.set_ylabel(ylabel, fontsize=fontsize_labels)
        elif col == 0:
            ax.set_ylabel("", fontsize=fontsize_labels)
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

        ax.xaxis.set_major_locator(plt.MaxNLocator(6))
        ax.yaxis.set_major_locator(plt.FixedLocator(y_ticks))
        
        legend = ax.legend(ncol=2, columnspacing=1, loc='upper left', labelspacing=0.4, fontsize=fontsize_textbox, frameon=False, handlelength=1, handletextpad=1)
        
        for handle in legend.legend_handles:
            handle.set_linewidth(4.0)
            
        ax.grid(False)
        ax.set_ylim(lower_limit, upper_limit)
        ax.set_xlim(-1, x_upper_limit)  # Setze die x-Achsen-Obergrenze
        ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)

    for j in range(len(ligands), num_rows * num_cols):
        fig.delaxes(axes.flatten()[j])

    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    plt.savefig(output_file, dpi=600)
    plt.close()


# Funktion zum Sammeln von .xvg-Dateien
import os

def collect_files(root_dir, proteins):
    file_dict = {protein: {'p1_rmsf': [], 'p2_rmsf': [], 'rmsf': [], 
                           'p1_rmsd': [], 'p2_rmsd': [], 'rmsd': [], 
                           'p1_gyration': [], 'p2_gyration': [], 'gyration': []} 
                 for protein in proteins}

    for root, _, files in os.walk(root_dir):
        for file in files:
            if file.endswith('.xvg'):
                # Beispiel: ./6FXb/CA_/6fxb_ca_1/prod_center_rottrans_p1_rmsf.xvg
                file_path = os.path.join(root, file)
                protein = next((p for p in proteins if p in root), None)
                
                if protein:
                    # Check file category and assign to the correct list
                    if 'rmsf' in file:
                        if 'p1' in file:
                            file_dict[protein]['p1_rmsf'].append(file_path)
                        elif 'p2' in file:
                            file_dict[protein]['p2_rmsf'].append(file_path)
                        else:
                            file_dict[protein]['rmsf'].append(file_path)
                    elif 'rmsd' in file:
                        if 'p1' in file:
                            file_dict[protein]['p1_rmsd'].append(file_path)
                        elif 'p2' in file:
                            file_dict[protein]['p2_rmsd'].append(file_path)
                        else:
                            file_dict[protein]['rmsd'].append(file_path)
                    elif 'gyration' in file:
                        if 'p1' in file:
                            file_dict[protein]['p1_gyration'].append(file_path)
                        elif 'p2' in file:
                            file_dict[protein]['p2_gyration'].append(file_path)
                        else:
                            file_dict[protein]['gyration'].append(file_path)
    
    return file_dict

# Funktion zur Erstellung individueller Plots
def plot_individual(protein, file_dict, output_prefix, colors=None):
    pairs = [
        ('p1_rmsf', 'p2_rmsf', 'RMSF', 'Aminosäure'),
        ('p1_rmsd', 'p2_rmsd', 'RMSD', 'Zeit [ns]'),
        ('p1_gyration', 'p2_gyration', 'Trägheitsradius', 'Zeit [ns]')
    ]
    
    for p1_category, p2_category, title, x_label in pairs:
        p1_files = file_dict[protein][p1_category]
        p2_files = file_dict[protein][p2_category]
        
        if p1_files or p2_files:
            plot_data(
                p1_files, 
                p2_files, 
                f'{title} von {protein.upper()}', 
                f'{output_prefix}_{title.upper()}.png', 
                ylabel=f'{title} [nm]', 
                x_label=x_label,
                colors=colors
            )


def read_xvg(file_path):
    x = []
    y_values = []
    # Bestimme, ob es sich um RMSF, RMSD oder Gyration handelt
    file_name_upper = file_path.upper()
    is_rmsf = 'RMSF' in file_name_upper
    is_rmsd = 'RMSD' in file_name_upper
    is_gyration = 'GYRATION' in file_name_upper
    
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('@') and not line.startswith('#'):
                data = line.split()
                if len(data) > 1:
                    # Wende die Zeitkonversion nur für RMSD und Gyration an
                    if not is_rmsf:  # Für RMSF keine Division
                        x.append(float(data[0]) / 1000.0)  # Convert ps to ns
                    else:
                        x.append(float(data[0]))  # Keine Umwandlung für RMSF

                    if not y_values:
                        y_values = [[] for _ in range(len(data) - 1)]
                    for i in range(len(data) - 1):
                        y_values[i].append(float(data[i + 1]))

    return np.array(x), np.array(y_values).T


def main():
    root_dir = '.'  # Ersetze dies durch dein Wurzelverzeichnis, falls anders
    proteins = ['6fxb']


    file_dict = collect_files(root_dir, proteins)


    for protein in proteins:
        plot_individual(protein, file_dict, protein.upper())

if __name__ == "__main__":
    main()
