#!/bin/bash

{
  echo "a1-2596"
  echo "name 24 Protein_p1"
  echo "a2597-5193"
  echo "name 25 Protein_p2"
  echo "24 & 2"
  echo "name 26 Protein-H_p1"
  echo "25 & 2"
  echo "name 27 Protein-H_p2"
  echo "1 | 13"
  echo "name 28 Protein_MOL"
  echo "q"
} | gmx make_ndx -f gromacs.gro -o index.ndx

echo 24| gmx genrestr -f gromacs.gro -n index.ndx -o posre_1000.itp -fc 1000
echo 24| gmx genrestr -f gromacs.gro -n index.ndx -o posre_500.itp -fc 500
echo 24| gmx genrestr -f gromacs.gro -n index.ndx -o posre_100.itp -fc 100
echo 24| gmx genrestr -f gromacs.gro -n index.ndx -o posre_10.itp -fc 10


