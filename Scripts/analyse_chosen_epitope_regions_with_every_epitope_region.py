import re
import os
from collections import defaultdict

alle_ueberschneidungen = []

lineare_epitope = """
ER4A, 1-16
ER4B, 31-60
ER4C, 67-86
ER4D, 127-152
ER5A, 58-77
ER5B, 76-95
ER5C, 121-140
ER6A, 1-15
ER6B, 56-72
ER6C, 76-90
ER6D, 136-150
ER6E, 6-20
ER6F, 111-125
ER7A, 41-60
ER7B, 102-124
ER7C, 149-162
"""

struktur_epitope = """
ER1A, Val3, Thr6, Ala80, Val81, Lys91, Val92, Leu93, Val94, Leu95, Leu104, Phe105, Cys106, Met107, Glu108
ER1B, Ile12, Gln13, Lys14, Val15, Ala16, Gly17, Thr18, Trp19, Tyr20, Val43, Glu44, Glu45, Leu46, Lys47, Pro48, Thr49, Pro50, Gly52, Asp53, Leu54, Glu55, Ile56, Leu57, Leu58, Lys70, Leu103, Leu122, Val123, Arg124		
ER1C, Leu22, Ala23, Met24, Ala25, Ala26, Leu32, Thr125, Pro126, Glu127, Val128, Asp129, Asp130, Leu133, Phe136, Asp137, Leu140, His146, Ile147, Arg148, Leu149, Ser150, Phe151, Asn152, Pro153, Leu156		
ER1C‘, Ser21, Leu22, Ala23, Met24, Gln35, Arg40, Tyr42, Val123, Thr125, Pro126, Glu127, Val128, Asp129, Leu149, Ser150, Phe151, Asn152, Pro153, Gln155, Leu156, Cys160, His161		
ER1D, Ser21, Leu22, Ala23, Met24, Ala25, Ala26, Ser27, Asp28, Ile29, Leu32, Arg40, Tyr42, Pro144, Met145, His146, Ile147, Arg148, Leu149, Ser150, Phe151, Asn152, Gln155, His161		
ER1E, Leu22, Lys101, Val123, Thr125, Pro126, Glu127, Val128, Asp129, Asp130, Glu131, Ala132, Leu133, Glu134, Lys135, Phe136, Asp137, Lys138, Pro153, Leu156
ER1F, Ile12, Lys47, Pro48, Thr49, Pro50, Glu51, Gly52, Asp53, Leu54, Glu55, Ile72, Ala73, Glu74, Lys75, Thr76, Val81, Phe82, Lys83, Ile84, Asp85		
ER1G, Ser27, Glu89, Asn90, Glu108, Asn109, Ser110, Ala111, Glu112, Pro113, Glu114, Gln115, Ser116, Leu117, Pro144		
ER1H, Gln5, Ala25, Leu95, Cys106, Cys119, Lys135, Phe136, Asp137, Lys138, Ala139, Leu140, Lys141, Ala142, Leu143, Pro144, Met145, His146, Arg148		
ER1I, Ala26, Ser27, Asp28, Ile29, Ser30, Leu31, Leu32, Asp33, Glu114, Gln115, Leu117, Pro144, Met145, His146, Ile147		
ER1J, Ile84, Asp85, Ala86, Leu87, Asn88, Glu89, Asn90, Asn109, Ala111, Glu112		
*ER2A, Thr18, Trp19, Tyr20, Val43, Glu44, Glu45, Leu46, Lys47, Leu57, Leu58, Gln59, Cys66, Ala67, Gln68, Pro126, Glu126, Thr154, Gln155, Leu156, Glu157
*ER3A, Ala139, Thr4, Val3, Ala142, Pro144, Ser27, Met145, Arg148, Leu149, Ser14
*ER3B, Leu117, Ser27, Met145, Pro144, Ala142, Ile2, Val3, Thr4, Gln5, Thr6, Met7, Lys8     
*ER3C, Leu46, Ala16, Pro48, Thr49, Gln13, Asp11, Lys14, Arg124, Thr125
*ER3D, Thr76, Lys77, Ile78, Pro79, Leu10, Asp11, Tyr99, Lys14, Gln13       
*ER3E, Thr4, Thr6, Met7, Lys8, Gly9, Asp11, Gln13, Lys14, Ala16, Arg124, Thr18, Tyr20      
*ER3F, Asn152, Thr154, Glu158, Glu157, Pro126, Thr125, Arg124, Tyr99, Leu10, Asp11, Gly9       
*ER3G, Glu158, Gln59, Gln159, Cys160, Cys66, Trp61, Gln35, Ser36, Ala37, Asp33, Ala34
*ER3H, Thr76, Lys75, Gly52, Thr49, Gln13, Ile12, Leu10, Asp11, Gly9, Lys8, Thr6
"""
#*Epitope has been -1 to position because the author had a systematic error

stable_contact_strings = """
2blg_clab_1, 1,9,12,49-53,72,74-79,81-89,91,110-111)
2blg_clab_2, 7,11-20,22,45-48,96-102,123-135,150-158)
2blg_clab_3, 4-9,14,22-24,29,32-35,40,42,95-96,98-102,124-135,137-139,141-142,147-158,160-162)
2blg_cla_1, 1-20,44-52,55,57,70,72,76-81,91,93,97-101,123-128,153-154,157)
2blg_cla_2, 1-14,48-53,74-84,89-96,99,106,108-113,138-139,141-144)
2blg_cad_1, 1-2,4-14,50-52,75-80,95-102,104,110,130-132,134-136,138-139,141-142)
2blg_ca_2_prod_mindistres_mol1_55_200_21200, 24,28-31,36-43,56,58-60,62,68-73,84-92,105-110,112-120,146)
2blg_clad_1, 1,7-14,48-54,72-91,93,99,108,110)
2blg_clad_2, 1-11,14,47-51,53,55,57,69-70,72,78-80,94-102,104,129-146,148)
6fxb_clad_1, 20,34-36,40,42,44,59,61,63-68,126-128,149-164,188-196,198-200,250,274-279,306-308)
6fxb_clad_2, 18-20,35,44-45,47,57-70,126,156-162)
6fxb_clad_3, 17-18,20,42-45,47-50,55-63,65-72,86,159-160,162)
6fxb_clab_1, 130,134,163-172,236-246,252-258,270-275,300,303-305)
6fxb_clab_2, 9-20,44-57,70-75,79,82-85,98-101,123-128,157)
6fxb_clab_3, 45,47-58,60,62-63,67-76,82-92,109)
6fxb_cla_1, 63-64,174-175,178-180,198-204,206-237,245-252,269,271,321)
6fxb_cla_3, 1-6,51,53,74-78,81-85,90-93,108-111)
6fxb_cad_1, 163-168,170,215,236-240,243-248,251-255,270-273,275)
6fxb_cad_2, 134,138,141,163-172,239-243,256-260,264,266-268,281,293,296-306)
6fxb_cad_3, 163-174,212,236-245,252-258,270,272,300,304)
6fxb_cab_1, 1-10,74-83,90-96,108,110,138)
6fxb_cab_2, 1-10,74-84,90-95,108-110,113)
6fxb_cab_3, 1-10,76-83,90-95,108-110)
"""

stable_contact_strings_2blg = """
2blg_clab_1, 1,9,12,49-53,72,74-79,81-89,91,110-111)
2blg_clab_2, 7,11-20,22,45-48,96-102,123-135,150-158)
2blg_clab_3, 4-9,14,22-24,29,32-35,40,42,95-96,98-102,124-135,137-139,141-142,147-158,160-162)
2blg_cla_1, 1-20,44-52,55,57,70,72,76-81,91,93,97-101,123-128,153-154,157)
2blg_cla_2, 1-14,48-53,74-84,89-96,99,106,108-113,138-139,141-144)
2blg_cad_1, 1-2,4-14,50-52,75-80,95-102,104,110,130-132,134-136,138-139,141-142)
2blg_ca_2_prod_mindistres_mol1_55_200_21200, 24,28-31,36-43,56,58-60,62,68-73,84-92,105-110,112-120,146)
2blg_clad_1, 1,7-14,48-54,72-91,93,99,108,110)
2blg_clad_2, 1-11,14,47-51,53,55,57,69-70,72,78-80,94-102,104,129-146,148)
"""

stable_contact_strings_6fxb = """
6fxb_clad_1, 20,34-36,40,42,44,59,61,63-68,126-128,149-164,188-196,198-200,250,274-279,306-308)
6fxb_clad_2, 18-20,35,44-45,47,57-70,126,156-162)
6fxb_clad_3, 17-18,20,42-45,47-50,55-63,65-72,86,159-160,162)
6fxb_clab_1, 130,134,163-172,236-246,252-258,270-275,300,303-305)
6fxb_clab_2, 9-20,44-57,70-75,79,82-85,98-101,123-128,157)
6fxb_clab_3, 45,47-58,60,62-63,67-76,82-92,109)
6fxb_cla_1, 63-64,174-175,178-180,198-204,206-237,245-252,269,271,321)
6fxb_cla_3, 1-6,51,53,74-78,81-85,90-93,108-111)
6fxb_cad_1, 163-168,170,215,236-240,243-248,251-255,270-273,275)
6fxb_cad_2, 134,138,141,163-172,239-243,256-260,264,266-268,281,293,296-306)
6fxb_cad_3, 163-174,212,236-245,252-258,270,272,300,304)
6fxb_cab_1, 1-10,74-83,90-96,108,110,138)
6fxb_cab_2, 1-10,74-84,90-95,108-110,113)
6fxb_cab_3, 1-10,76-83,90-95,108-110)
"""

einbuchstaben_zu_dreibuchstaben = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
}
blg_fasta = ("LIVTQTMKGLDIQKVAGTWYSLAMAASDISLLDAQSAPLRVYVEELKPTPEGDLEILLQKWENDECAQKKIIAEKTKIPAVFKIDALNENKVLVLDTDYKKYLLFCMENSAEPEQSLVCQCLVRTPEVDDEALEKFDKALKALPMHIRLSFNPTQLEEQCHI")


def generate_pymol_selection_for_most_common_amino_acids(amino_acid_positions, threshold):
    pymol_selections = []

    # Iteriere über die Aminosäuren und deren Positionen
    for amino_acid, positions in amino_acid_positions.items():
        for position, count in positions.items():
            if count >= threshold:
                pymol_selections.append(f"{position}")
    
    return pymol_selections

def generiere_pymol_selektion(liste):
    # Erster Eintrag ist der Name, Rest sind Aminosäuren
    name = liste[0]
    aminosäuren = liste[1:]
    
    # Extrahiere die Nummern aus den Aminosäuren
    nummern = []
    for eintrag in aminosäuren:
        # Suche die Zahl am Ende des Eintrags
        zahl = ''.join(filter(str.isdigit, eintrag))
        if zahl:  # nur hinzufügen, wenn eine Zahl gefunden wurde
            nummern.append(zahl)

    # Erstelle den PyMOL-Selektionsbefehl
    pymol_selektion = f"select {name}, res {'+'.join(nummern)}"
    return pymol_selektion


def umwandeln_und_nummerieren(sequence):
    # Erzeugt eine nummerierte dreibuchstabige Aminosäurensequenz
    dreibuchstaben_sequence = []
    for position, buchstabe in enumerate(sequence, start=1):
        dreibuchstaben_code = einbuchstaben_zu_dreibuchstaben.get(buchstabe, '')
        dreibuchstaben_sequence.append(f"{dreibuchstaben_code}{position}")
    return dreibuchstaben_sequence

def bereich_zu_aminosaeuren(substrings, fasta_sequenz):
    ergebnis = []
    
    for item in substrings:
        # Entferne alle nicht-numerischen oder nicht-Bereich-Zeichen
        item = item.strip().replace("res", "").replace(")", "").replace("(", "")
        
        if "-" in item:  # Bereich
            start, ende = map(int, item.split('-'))
            for pos in range(start, ende + 1):
                if pos > 162: pos -= 162
                aminosaeure = einbuchstaben_zu_dreibuchstaben.get(fasta_sequenz[pos - 1], "Unk")  # Verwende "Unk" für unbekannte Aminosäuren
                ergebnis.append(f"{aminosaeure}{pos}")
        
        elif item.isdigit():  # Einzelne Position (eine Zahl ohne Aminosäurenamen)
            pos = int(item)
            if pos > 162: pos -= 162
            aminosaeure = einbuchstaben_zu_dreibuchstaben.get(fasta_sequenz[pos - 1], "Unk")
            ergebnis.append(f"{aminosaeure}{pos}")
        
        else:  # Bereits ein Aminosäurencode (z.B., Thr77)
            ergebnis.append(item.strip())

    return ergebnis

def count_amino_acids_with_positions(strings_list):

    amino_acid_positions = defaultdict(lambda: defaultdict(int))
    
    # Regex-Muster, um Aminosäuren mit Positionen zu extrahieren
    amino_acid_pattern = r"\b([A-Z][a-z]{2})(\d+)\b"

    # Verarbeitung der Liste
    for line in strings_list:
        # Überspringe `None`-Werte und nicht-String-Einträge
        if line is None or not isinstance(line, str):
            continue
        
        # Finde alle Matches (Aminosäure + Position)
        matches = re.findall(amino_acid_pattern, line)
        
        # Zähle die Vorkommen der Positionen pro Aminosäure
        for amino_acid, position in matches:
            amino_acid_positions[amino_acid][int(position)] += 1
    
    return amino_acid_positions



# Funktion zum Parsen und Umwandeln
def extrahiere_aminosaeuren_aus_bereich(bereich_string):
    substring = [teil.strip(' ') for teil in bereich_string.split(',')]
    return bereich_zu_aminosaeuren(substring, blg_fasta)



def sortiere_nach_residuen(liste):
    """Sortiert eine Liste von Aminosäuren anhand ihrer Residuen."""
    # Angenommen, jedes Element in der Liste hat die Form "AminosäureX" (z.B. 'Thr76')
    return sorted(liste, key=lambda x: int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else float('inf'))

def finde_ueberschneidungen(liste1, liste2):
    """Findet die Überschneidungen zwischen zwei Listen und gibt sie formatiert aus."""
    

    # Sortiere beide Listen basierend auf den Zahlen im Namen
    sortierte_liste1 = sortiere_nach_residuen(liste1[1:])
    sortierte_liste2 = sortiere_nach_residuen(liste2[1:])
    
    # Extrahiere den Namen der Listen
    name1 = liste1[0]
    name2 = liste2[0]
    
    # Finde die Überschneidungen
    ueberschneidungen = set(sortierte_liste1).intersection(sortierte_liste2)
    
    # Ausgabe formatieren
    if ueberschneidungen:
    #    print(name2, end=', ')
    #    print(f"{name1} überschneidet {name2}")
        #print(f"{name1} überschneidet {name2} mit: {', '.join(sorted(ueberschneidungen, key=lambda x: int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else float('inf')))}")
        return (f"{name1} überschneidet {name2} mit: {', '.join(sorted(ueberschneidungen, key=lambda x: int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else float('inf')))}")
    #else:
    #    print(f"Keine Überschneidungen zwischen {name1} und {name2}")



print("\n", "VERGLEICHE STABILE BLG-LIGANDEN-KONTAKTE MIT ALLEN EPITOPREGIONEN:\n")

#ph6 
# print alle überschneidungen der stabilen liganden und bekannten epitope. printe das ergebnis und speichere das ergebnis in string cache
string_cache = []
for line1 in stable_contact_strings_6fxb.splitlines():
    print("")
    #print(extrahiere_aminosaeuren_aus_bereich(line1)[0])
    extrahiere_aminosaeuren_aus_bereich(line1)[0]
    for line2 in struktur_epitope.splitlines():
        string_cache.append(finde_ueberschneidungen(extrahiere_aminosaeuren_aus_bereich(line1),extrahiere_aminosaeuren_aus_bereich(line2)))
    for line2 in lineare_epitope.splitlines():
        string_cache.append(finde_ueberschneidungen(extrahiere_aminosaeuren_aus_bereich(line1),extrahiere_aminosaeuren_aus_bereich(line2)))

# parse durch string cache und zähle die aminosäuren und sortiere sie nach häufigkeit
amino_acid_positions = count_amino_acids_with_positions(string_cache)
sorted_amino_acids = []
for amino_acid, positions in amino_acid_positions.items():
    for position, count in positions.items():
        sorted_amino_acids.append((amino_acid, position, count))
sorted_amino_acids.sort(key=lambda x: (-x[2], x[0], x[1]))
for amino_acid, position, count in sorted_amino_acids:
    print(f"{amino_acid}[{position}]: {count}", end = ', ')

# setze den threshold der häufigkeit für die pymol selektion und generiere die pymol selektions ausgabe
threshold=0
pymol_selections = generate_pymol_selection_for_most_common_amino_acids(amino_acid_positions, threshold)
print("PyMOL Selektionskommandos 6fxb: select most_common_AA_6fxb_t0, res ", end='')
for selection in pymol_selections:
    print(selection, end='+')
threshold=12
pymol_selections = generate_pymol_selection_for_most_common_amino_acids(amino_acid_positions, threshold)
print("PyMOL Selektionskommandos 6fxb: select most_common_AA_6fxb_t12, res ", end='')
for selection in pymol_selections:
    print(selection, end='+')
threshold=20
pymol_selections = generate_pymol_selection_for_most_common_amino_acids(amino_acid_positions, threshold)
print("PyMOL Selektionskommandos 6fxb: select most_common_AA_6fxb_t20, res ", end='')
for selection in pymol_selections:
    print(selection, end='+')

#ph9 
# print alle überschneidungen der stabilen liganden und bekannten epitope. printe das ergebnis und speichere das ergebnis in string cache
string_cache = []
for line1 in stable_contact_strings_2blg.splitlines():
    print("")
    #print(extrahiere_aminosaeuren_aus_bereich(line1)[0])
    extrahiere_aminosaeuren_aus_bereich(line1)[0]
    for line2 in struktur_epitope.splitlines():
        string_cache.append(finde_ueberschneidungen(extrahiere_aminosaeuren_aus_bereich(line1),extrahiere_aminosaeuren_aus_bereich(line2)))
    for line2 in lineare_epitope.splitlines():
        string_cache.append(finde_ueberschneidungen(extrahiere_aminosaeuren_aus_bereich(line1),extrahiere_aminosaeuren_aus_bereich(line2)))

# parse durch string cache und zähle die aminosäuren und sortiere sie nach häufigkeit
amino_acid_positions = count_amino_acids_with_positions(string_cache)
sorted_amino_acids = []
for amino_acid, positions in amino_acid_positions.items():
    for position, count in positions.items():
        sorted_amino_acids.append((amino_acid, position, count))
sorted_amino_acids.sort(key=lambda x: (-x[2], x[0], x[1]))
for amino_acid, position, count in sorted_amino_acids:
    print(f"{amino_acid}[{position}]: {count}", end = ', ')

# setze den threshold der häufigkeit für die pymol selektion und generiere die pymol selektions ausgabe
threshold=0
pymol_selections = generate_pymol_selection_for_most_common_amino_acids(amino_acid_positions, threshold)
print("PyMOL Selektionskommandos 2blg: select most_common_AA_2blg_t0, res ", end='')
for selection in pymol_selections:
    print(selection, end='+')
# setze den threshold der häufigkeit für die pymol selektion und generiere die pymol selektions ausgabe
threshold=12
pymol_selections = generate_pymol_selection_for_most_common_amino_acids(amino_acid_positions, threshold)
print("PyMOL Selektionskommandos 2blg: select most_common_AA_2blg_t12, res ", end='')
for selection in pymol_selections:
    print(selection, end='+')
# setze den threshold der häufigkeit für die pymol selektion und generiere die pymol selektions ausgabe
threshold=20
pymol_selections = generate_pymol_selection_for_most_common_amino_acids(amino_acid_positions, threshold)
print("PyMOL Selektionskommandos 2blg: select most_common_AA_2blg_t20, res ", end='')
for selection in pymol_selections:
    print(selection, end='+')

print("\n", "VERGLEICHE EPITOPREGIONEN 1, 2 UND 3 MIT ALLEN EPITOPREGIONEN:\n")


"""


for line2 in struktur_epitope.splitlines():
    print("")
    finde_ueberschneidungen(extrahiere_aminosaeuren_aus_bereich('ER7A, 41-60'),extrahiere_aminosaeuren_aus_bereich(line2))
for line2 in lineare_epitope.splitlines():
    print("")
    finde_ueberschneidungen(extrahiere_aminosaeuren_aus_bereich('ER7A, 41-60'),extrahiere_aminosaeuren_aus_bereich(line2))

for line2 in struktur_epitope.splitlines():
    print("")
    finde_ueberschneidungen(extrahiere_aminosaeuren_aus_bereich('ER7B, 102-124'),extrahiere_aminosaeuren_aus_bereich(line2))
for line2 in lineare_epitope.splitlines():
    print("")
    finde_ueberschneidungen(extrahiere_aminosaeuren_aus_bereich('ER7B, 102-124'),extrahiere_aminosaeuren_aus_bereich(line2))
    
for line2 in lineare_epitope.splitlines():
    print("")
    finde_ueberschneidungen(extrahiere_aminosaeuren_aus_bereich('ER7C, 149-162'),extrahiere_aminosaeuren_aus_bereich(line2))       
for line2 in struktur_epitope.splitlines():
    print("")
    finde_ueberschneidungen(extrahiere_aminosaeuren_aus_bereich('ER7C, 149-162'),extrahiere_aminosaeuren_aus_bereich(line2))



for line2 in struktur_epitope.splitlines():
    print(generiere_pymol_selektion(extrahiere_aminosaeuren_aus_bereich(line2)))

print(extrahiere_aminosaeuren_aus_bereich("41-60"))
print(extrahiere_aminosaeuren_aus_bereich("102-124"))
print(extrahiere_aminosaeuren_aus_bereich("149-162"))
"""
print(generiere_pymol_selektion(extrahiere_aminosaeuren_aus_bereich("ER7A, 41-60")))
print(generiere_pymol_selektion(extrahiere_aminosaeuren_aus_bereich("ER7B, 102-124")))
print(generiere_pymol_selektion(extrahiere_aminosaeuren_aus_bereich("ER7C, 149-162")))


