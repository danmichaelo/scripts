#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

from ase.io.vasp import read_vasp, write_vasp
atoms = read_vasp('POSCAR')
transvec = atoms.get_cell().diagonal() / -1.
atoms.translate(transvec)
supercell = atoms.repeat((3,3,3))
write_vasp('POSCAR.3',supercell)


 
# import numpy as np
# from oppvasp.vasp.parsers import PoscarParser
# 
# pp = PoscarParser('POSCAR')
# unit_pos = pp.get_positions( coordinates = 'direct' )
# unit_natoms = pos.shape[0]
# print "Num atoms:",num_atoms
# super_size = [2,2,2]
# natoms = np.prod(super_size)*num_atoms
# pos = np.array((natoms,3))
# pos[0:unit_natoms] = unit_pos
# 
# # add atoms in x dir:
# for i = range(unit_natoms):
#     for j = range(1,super_size[0]):
#         pos[unit_natoms+i*j] = pos[i]*[j,1,1]

# f = open('movie.xyz','w')
# for i in range(ndsteps):
#     f.write("%d\n" % (nions)) 
#     f.write("Energy: %.3f\n" % data['e_fr_energy'][i])
#     basis = data['basis'][i]
#     for j in range(nions):
#         # convert to cartesian (Angstrom):
#         at_pos = np.dot(data['atoms'][j]['trajectory'][i], basis)
#         f.write(atoms[j] + '    ' + "%.8f    %.8f    %.8f\n" % tuple(at_pos))
# f.close()
# 
# # Print final position of the first two atoms in direct coordinates:
# print data['atoms'][0]['trajectory'][ndsteps-1]
# print data['atoms'][1]['trajectory'][ndsteps-1]
# 
