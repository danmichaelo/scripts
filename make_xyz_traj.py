#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import numpy as np
from oppvasp.vasp.parsers import IterativeVasprunParser

vp = IterativeVasprunParser('vasprun.xml')
nsteps = vp.get_num_ionic_steps() # NSW
nions = vp.get_num_atoms()
atoms = vp.get_atoms()
data = vp.get_all_trajectories()
ndsteps = data['length'] # number of steps actually calculated

f = open('movie.xyz','w')
for i in range(ndsteps):
    f.write("%d\n" % (nions)) 
    f.write("Energy: %.3f\n" % data['e_fr_energy'][i])
    basis = data['basis'][i]
    for j in range(nions):
        # convert to cartesian (Angstrom):
        at_pos = np.dot(data['atoms'][j]['trajectory'][i], basis)
        f.write(atoms[j] + '    ' + "%.8f    %.8f    %.8f\n" % tuple(at_pos))
f.close()

# Print final position of the first two atoms in direct coordinates:
print data['atoms'][0]['trajectory'][ndsteps-1]
print data['atoms'][1]['trajectory'][ndsteps-1]

