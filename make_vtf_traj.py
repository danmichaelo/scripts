#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import os,sys
import numpy as np
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
             RotatingMarker, ReverseBar, SimpleProgress

from oppvasp import getAtomicNumberFromSymbol
from oppvasp.vasp.parsers import IterativeVasprunParser, PoscarParser
from oppvasp.md import Trajectory

xml_filename = 'vasprun.xml'
npz_filename = 'trajectory.npz'
vtf_filename = 'trajectory.vtf'

if len(sys.argv) > 1:
    npz_filename = sys.argv[1]

if os.path.isfile(npz_filename):
    traj = Trajectory(filename = npz_filename)
else:
    p = IterativeVasprunParser(xml_filename)
    traj = p.get_all_trajectories()
    traj.save(npz_filename)

basis = traj.basis
atoms = traj.atoms
nsteps = traj.length # NSW
nions = traj.num_atoms
positions = traj.positions
forces = traj.forces
total_energy = traj.total_energy

f = open(vtf_filename,'w')

# A bond is formed whenever two atoms are within (R1 + R2) x 0.6 of each other, 
# where R1 and R2 are the respective radii of candidate atoms.
R = 2.0

def writeAtomDef(b,e,c):
    if b == e-1:
        str = 'atom %d' % (b)
    else:
        str = 'atom %d:%d' % (b,e)
    str += ' radius %.1f atomicnumber %d\n' % (R,c)
    f.write(str)


f.write("# Structure block\n")
b=0; e=0; c=-1
for at in atoms:
    if c != at and e > 0:
        writeAtomDef(b,e,c)
        b = e
    c = at
    e += 1
writeAtomDef(b,e-1,c)

def norm(v):
    return np.sqrt(np.dot(v,v))

def getUnitCell(i):
    A=basis[i][0]
    B=basis[i][1]
    C=basis[i][2]
    a = norm(A)
    b = norm(B)
    c = norm(C)
    alpha = np.arccos( np.dot(B,C) / (b*c) ) * 360 / (2.*np.pi)
    beta  = np.arccos( np.dot(C,A) / (c*a) ) * 360 / (2.*np.pi)
    gamma = np.arccos( np.dot(A,B) / (a*b) ) * 360 / (2.*np.pi)
    return a,b,c,alpha,beta,gamma

pbar = ProgressBar(widgets=[u'الرجاء الانتظار',Bar()], maxval = nsteps).start()
for i in range(nsteps):
    pbar.update(i)
    f.write("\n# Energy: %.3f\n" % total_energy[i])
    f.write("timestep\n")
    uc = getUnitCell(i)
    f.write("unitcell %.2f %.2f %.2f %.1f %.1f %.1f\n" % uc)
    for j in range(nions):
        # convert to cartesian (Angstrom):
        at_pos = np.dot(positions[i,j], basis[i])
        at_force = forces[i,j]
        # should probably looking into optimization here
        f.write("%.8f    %.8f    %.8f     %.8f    %.8f    %.8f\n" % tuple(list(at_pos) + list(at_force)))
pbar.finish()
f.close()
print "Wrote %s" % (vtf_filename)
