#!/bin/bash

echo 2| gmx genrestr -f gromacs.gro -o posre_1000.itp -fc 1000
echo 2| gmx genrestr -f gromacs.gro -o posre_500.itp -fc 500
echo 2| gmx genrestr -f gromacs.gro -o posre_100.itp -fc 100
echo 2| gmx genrestr -f gromacs.gro -o posre_10.itp -fc 10

{
  echo "1 | 13"
  echo "name 24 Protein_MOL"
  echo "2 | 13"
  echo "name 25 Protein-H_MOL"
  echo "q"
} | gmx make_ndx -f gromacs.gro -o index.ndx
