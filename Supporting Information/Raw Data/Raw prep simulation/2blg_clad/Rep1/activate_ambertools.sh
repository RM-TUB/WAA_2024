
#!/bin/bash

cp "./../EFZ.frcmod" "./EFZ.frcmod"
cp "./../EFZ_sybyl_resp_gaff.prepc" "./EFZ_sybyl_resp_gaff.prepc"

# Aktiviere die AmberTools23-Umgebung
source activate AmberTools23

# Führe tleap aus
tleap -f tleap.in
