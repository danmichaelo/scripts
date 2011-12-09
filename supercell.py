#!/usr/bin/env python

import sys, os
import numpy as np
from oppvasp.vasp.parsers import PoscarParser
from oppvasp.structure import Structure 
from oppvasp.utils import query_yes_no

## Parameters
infile = 'CONTCAR.1'
outfile = 'CONTCAR-222'
px = [-1,0]
py = [-1,0]
pz = [-1,0]


try:
    pp = PoscarParser(infile)
except:
    print "Failed to read %s" % (infile)
    sys.exit(1)
s = pp.get_structure()
p = s.get_positions('d')
a = s.get_atomtypes()
c = s.get_cell()

nx = len(px)
ny = len(py)
nz = len(pz)

singlesize = p.shape[0]
totalsize = singlesize * len(px)*len(py)*len(pz)
print "Original unit cell contained %d atoms. New cell will contain %d atoms" % (singlesize,totalsize)

p2 = np.zeros((totalsize,3))
a2 = np.zeros((totalsize))

for i,x in enumerate(px):
    for j,y in enumerate(py):
        for k,z in enumerate(pz):
            f = len(py)*len(pz)*i + len(pz)*j + k
            print "Adding cell %d of %d" % (f+1,nx*ny*nz)
            ids = f*singlesize
            ide = (f+1)*singlesize
            p2[ids:ide,0] = (x + p[:,0])  / nx
            p2[ids:ide,1] = (y + p[:,1])  / ny
            p2[ids:ide,2] = (z + p[:,2])  / nz
            a2[ids:ide] = a

c2 = np.zeros(c.shape)
c2[0] = c[0] * nx
c2[1] = c[1] * ny
c2[2] = c[2] * nz

s2 = Structure( cell = c2, positions = p2, atomtypes = a2, coords = 'direct')
if os.path.exists(outfile):
    print "\nWARNING: The file \"%s\" already exists." % outfile
    if query_yes_no("         Do you want to overwrite it?") == 'no':
        print "Ok, no supercell saved. Exiting..."
        sys.exit(0)
s2.save(outfile)
print "Saved %dx%dx%d cell to %s" % (nx,ny,nz,outfile)

