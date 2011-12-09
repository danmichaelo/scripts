#!/usr/bin/env python
# encoding=utf8

import sys
import numpy as np
from ase.io import vasp

# Assuming cubic cell !!!

filename = 'POSCAR'
if len(sys.argv) > 1:
    filename = sys.argv[1]

print "Checking",filename

v = vasp.read_vasp(filename)
r0 = v.get_scaled_positions()
a = v.get_cell()[0,0]

natoms = r0.shape[0]
print "Found",natoms,"atoms, lattice parameter is",a

# Shift all positions
x = r0[:,0]
minx = np.min(x[x>0])
invx = 1./minx
print "minx is",minx,"inverse is",invx
kvec = np.array((invx,0,0)) # invers av 1/(4s) der s er supercellestørrelsen

kvec = np.array((8,0,0)) # invers av 1/(4s) der s er supercellestørrelsen

smallshift = 1./a*0.1 # 0.2 Å
largeshift = 1./a*2. # 2 Å

theta = np.pi/5
rotmat = np.array((
    (np.cos(theta), np.sin(theta), 0),
    (-np.sin(theta), np.cos(theta), 0),
    (0, 0, 1)
    ))
R = [
    ('original', r0),
    ('uniform shift of all atoms', r0 + 1/(4*invx)),
    ('uniform rotation of all atoms', np.dot(r0, rotmat)),
    ('small random shift (± 0.2 Å) of all atoms', r0 + (-smallshift+np.random.rand(natoms,3)*2*smallshift)),
    ('large random shift (± 2 Å)', r0 + (-largeshift+np.random.rand(natoms,3)*2*largeshift))
]


for desc, r in R:
    print "Testing:",desc

    #move inside PBCs
    r -= np.floor(r)

    fac = 2 * np.pi * kvec[0] 
    print " -> T1(x):", np.sum(np.cos(fac * r[:,0]))/natoms 
    print " -> T1:", 1./3 * (np.sum(np.cos(fac * r[:,0]))/natoms + np.sum(np.cos(fac * r[:,1]))/natoms + np.sum(np.cos(fac * r[:,2]))/natoms )

    vd = 2*np.pi*np.dot(r, kvec)

    print " -> T2:", np.sum(np.exp(np.complex(0,1) * vd)).real / natoms

    cossum = np.sum(np.cos(vd))
    sinsum = np.sum(np.sin(vd))
    cossumX = np.sum(np.cos(fac*r[:,0]))
    sinsumX = np.sum(np.sin(fac*r[:,0]))
    print " -> T3(x):", np.sqrt(cossumX**2+sinsumX**2)/natoms," (origin independent)"
    print " -> T3:", np.sqrt(cossum**2+sinsum**2)/natoms," (origin independent)"

print "In a disordered system, T3 should be of order sqrt(natoms) =",np.sqrt(1./natoms)


