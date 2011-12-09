#!/usr/bin/env python
#encoding=utf8

import sys,os
import numpy as np
import argparse 
from oppvasp.vasp.parsers import PoscarParser
import oppvasp.util 
from oppvasp import direct_to_cartesian

parser = argparse.ArgumentParser( description = 'Adds small random displacements to all atoms, and saves as a new file.' )

parser.add_argument('--max', '-m', type=float, default=0.4, help='Max displacement, in Angstrom' )
#parser.add_argument('--no_wrap', '-n', action='store_true', help='Don\'t wrap coordinates back into unit cell' )
parser.add_argument('infile', nargs='?', type=oppvasp.util.FileType('r'), default='POSCAR', help='Input filename (defaults to POSCAR)')
parser.add_argument('outfile', nargs='?', type=oppvasp.util.FileType('w'), default=sys.stdout, help='Output filename (defaults to stdout)')
args = parser.parse_args()

# Read atoms from POSCAR:
structure1 = PoscarParser(args.infile, silent=True).get_structure()
pos = structure1.get_positions( coords = 'cart' )

#print "Max displacement: %.3f Angstrom" % args.max
# find random displacements using spherical coordinates:
rnd = np.random.random(pos.shape)
r = rnd[:,0] * args.max
azimuth = rnd[:,1]*2*np.pi
zenith = rnd[:,1]*np.pi

# convert to cartesian coordinates:
rnd[:,0] = r * np.cos(azimuth) * np.sin(zenith)
rnd[:,1] = r * np.sin(azimuth) * np.sin(zenith)
rnd[:,2] = r * np.cos(zenith)

structure1.set_positions( pos + rnd , coords='cart' )

structure1.save(args.outfile, direct_coords = True)

