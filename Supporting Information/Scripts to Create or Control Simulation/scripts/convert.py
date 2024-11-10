import parmed as pmd
amber = pmd.load_file('ref.prmtop', 'ref.inpcrd')
amber.save('topol.top')
amber.save('gromacs.gro')

