#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import os,sys
import numpy as np

try:
    from vtk_xml_serial_unstructured import VTK_XML_Serial_Unstructured
except ImportError:
    print "Error: The Python module 'vtk_xml_serial_unstructured' is not available. "
    print "It can (hopefully) be downloaded from "
    print "http://www.shocksolution.com/wordpress/wp-content/uploads/2009/01/vtk_xml_serial_unstructured-01.py"
    sys.exit(1)

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
total_energy = traj.total_energy

vtk = VTK_XML_Serial_Unstructured()

# A bond is formed whenever two atoms are within (R1 + R2) x 0.6 of each other, 
# where R1 and R2 are the respective radii of candidate atoms.
R = 2.0

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

for i in range(nsteps):
    filename = 'vtk_%d.vtu' % (i)
    vtk.snapshot(filename, positions[i,:,0], positions[i,:,1], positions[i,:,2])

vtk.writePVD('vtk.pvd')




