#!/usr/bin/env python
#encoding=utf8

import sys
import numpy as np
import argparse 
from oppvasp.vasp.parsers import PoscarParser
from oppvasp.util import get_pairs
from oppvasp import direct_to_cartesian

parser = argparse.ArgumentParser( description = 'Compares two POSCAR-type files, and prints a summary' )

parser.add_argument('--diff', '-d', action='store_true',
        help='Print displacement vectors for all atoms' )
parser.add_argument('--abs_diff', '-a', action='store_true',
        help='Print abs(displacements) for all atoms' )
parser.add_argument('--rad_diff', '-r', action='store_true',
        help='Print radial displacements (displacement vector norms) for all atoms' )
parser.add_argument('--no_unwrap', '-n', action='store_true',
        help='Don\'t try to unwrap motion over periodic boundaries.' )
parser.add_argument('infile', nargs='?', default='POSCAR', type=argparse.FileType('r'), help='POSCAR filename')
args = parser.parse_args()

# Read atoms from POSCAR:
poscar1 = PoscarParser(args.infile).get_structure()
pos = poscar1.get_positions( coords = 'direct' )
natoms = pos.shape[0]

# Build array of pair indices:
pairs = get_pairs(natoms)
npairs = pairs.shape[0]

print "%d atoms, %d pairs" % (natoms, npairs)

# Find displacement vectors for all atoms:
x = pos[pairs[:,0]] - pos[pairs[:,1]]

# Use minimum image convention to threat bonds over PBCs
# Note: This will not work with *very* tilted unit cells
x = x - (2*x).astype('int')

X = direct_to_cartesian(x, poscar1.get_cell())

r2 = (X**2).sum(axis=1)
r = np.sqrt(r2)

meanr2 = np.mean(r2,axis=0)
meanr = np.mean(r,axis=0)

minidx = np.argmin(r)
minpair = pairs[minidx]
minr = r[minidx]

print "Shortest bond is: %.3f Angstrom (between atoms %d and %d)" % (minr,minpair[0],minpair[1])
