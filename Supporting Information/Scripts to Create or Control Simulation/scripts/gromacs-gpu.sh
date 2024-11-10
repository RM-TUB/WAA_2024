#!/bin/bash --login
#SBATCH --job-name=gpu
#SBATCH --output=gpu.%N.%A.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --time=124:00:00


## check mdp files - same temperature and pressure of experiments?

module rm ompi
module load gromacs/2020.1

## -maxwarn removed 


# energy minimization
gmx grompp -f em.mdp -c gromacs.gro -p topol.top -n index.ndx -o em.tpr -r gromacs.gro  # CHANGE name of file.gro 
gmx mdrun -deffnm em

# heating
gmx grompp -f heat1.mdp -c em.gro -p topol.top -n index.ndx -o heat1.tpr -r gromacs.gro
gmx mdrun -deffnm heat1 

let a=2
while (($a <= 4))
do
        let prev=$a-1
        gmx grompp -f heat"$a".mdp -c heat"$prev".gro -p topol.top -n index.ndx -o heat"$a".tpr -r gromacs.gro
        gmx mdrun -deffnm heat"$a" 
        let a=$a+1
done

# pressure adjustment
gmx grompp -f press.mdp -c heat4.gro -p topol.top -n index.ndx -o press.tpr -r gromacs.gro -maxwarn 1
gmx mdrun -deffnm press 

# removal of positional restraints
gmx grompp -f pos1.mdp -c press.gro -p topol.top -n index.ndx -o pos1.tpr -r gromacs.gro -maxwarn 1
gmx mdrun -deffnm pos1 

let a=2
while (($a <= 4))
do
        let prev=$a-1
        gmx grompp -f pos"$a".mdp -c pos"$prev".gro -p topol.top -n index.ndx -o pos"$a".tpr -r gromacs.gro -maxwarn 1
        gmx mdrun -deffnm pos"$a" 
        let a=$a+1
done

# productive simulation
gmx grompp -f prod.mdp -c pos4.gro -p topol.top -n index.ndx -o prod.tpr -r gromacs.gro 
gmx mdrun -deffnm prod 


