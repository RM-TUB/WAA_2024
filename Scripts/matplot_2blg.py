import os
import matplotlib.pyplot as plt
import numpy as np
import colorsys
import re
import matplotlib.gridspec as gridspec

def extract_replica_number(file_path):
    match = re.search(r'_cab_(\d+)', file_path)  # Beispielregex, um die Zahl nach "_ca_" zu extrahieren
    if match:
        return int(match.group(1))
    return 0

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
                    x.append(float(data[0]))  # ps zu ns
                    if not y_values:
                        y_values = [[] for _ in range(len(data) - 1)]
                    for i in range(len(data) - 1):
                        y_values[i].append(float(data[i + 1]))
    return np.array(x), np.array(y_values).T


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


def plot_data(files, title, output_file, data_index=0, ylabel='Value [nm]', colors=None, alpha=0.7, x_label='Time [ns]', legend_labels=None, custom_max_y=None, constrained_layout=True):

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

    global_y_min, global_y_max = float('inf'), float('-inf')
    for file_path in files:
        x, y_values = parse_xvg(file_path)
        global_y_min = min(global_y_min, min(y_values[:, data_index]))
        global_y_max = max(global_y_max, max(y_values[:, data_index]))

    if 'RMSF' in title.upper():
        upper_limit = 0.64
        lower_limit = 0
        x_upper_limit = 163

    elif 'RMSD' in title.upper():
        upper_limit = 0.36
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

        ligand_files = [f for f in files if ligand in f]

        ligand_files.sort(key=extract_replica_number)
        
        for j, file in enumerate(ligand_files):
            alpha = [1.0, 0.5, 0.25]
            x, y = read_xvg(file)
            graphs = ax.plot(x, y[:, data_index], color=colors[ligand], alpha=alpha[j], label=f'BLG-{legend_labels[ligand]} ({j+1})', linewidth=0.5)
            if 'RMSF' in title.upper():
                for graph in graphs:
                    graph.set_linewidth(1)
        ligand_files = [f for f in files if 'NOLIGAND' in f]
        
        ligand_files.sort(key=extract_replica_number)

        for idx, p1_file in enumerate(ligand_files):
            alpha = [0.6, 0.3, 0.15]
            x_p1, y_p1 = read_xvg(p1_file)
            label=f'BLG-NK ({idx+1})'
            graphs = ax.plot(x_p1, y_p1[:, data_index], color=colors['NOLIGAND'], alpha=alpha[idx], label=label, linewidth=0.5)
            if 'RMSF' in title.upper():
                for graph in graphs:
                    graph.set_linewidth(1)        
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

    plt.subplots_adjust(hspace=0.1, wspace=0.1)  # hspace für vertikale und wspace für horizontale Abstände
    plt.savefig(output_file, dpi=600)
    plt.close()

def plot_individual(protein, file_dict, output_prefix, colors=None, custom_max_y_rmsf=0.44):
    """Plot individual data for a protein in subplots."""
    categories = [
        ('rmsf', 'RMSF', 'Aminosäure'),
        ('rmsd', 'RMSD', 'Zeit [ns]'),
        ('gyration', 'Gyration', 'Zeit [s]')
    ]
    
    for category, title, x_label in categories:
        files = file_dict[protein][category]
        
        if files:
            plot_data(
                files, 
                f'{title} of {protein.upper()}', 
                f'{output_prefix}_{title.upper()}.png', 
                ylabel=f'{title} (nm)', 
                x_label=x_label,
                colors=colors,  
                custom_max_y=custom_max_y_rmsf if category == 'rmsf' else None
            )

def collect_files(root_dir, proteins):
    """Collect .xvg files from directories."""
    file_dict = {protein: {'rmsf': [], 'rmsd': [], 'gyration': []} for protein in proteins}
    
    for root, _, files in os.walk(root_dir):
        for file in files:
            if file.endswith('.xvg'):
                file_path = os.path.join(root, file)
                for protein in proteins:
                    if protein in root:
                        if 'prod_center_rottrans_rmsf.xvg' in file:
                            file_dict[protein]['rmsf'].append(file_path)
                        elif 'prod_center_rottrans_rmsd.xvg' in file:
                            file_dict[protein]['rmsd'].append(file_path)
                        elif 'prod_center_rottrans_gyration.xvg' in file:
                            file_dict[protein]['gyration'].append(file_path)
    return file_dict

def main():
    root_dir = '.'  # Replace with your root directory if different
    proteins = ['2blg']

    file_dict = collect_files(root_dir, proteins)

    for protein in proteins:
        plot_individual(protein, file_dict, protein.upper())

if __name__ == "__main__":
    main()
