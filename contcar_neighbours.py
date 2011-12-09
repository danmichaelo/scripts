#!/usr/bin/env python
# encoding=utf8

import numpy as np
import ase
import sys
from ase.calculators.neighborlist import NeighborList

class BondList:

    def __init__(self, filename):
        atoms = ase.io.read(filename, format='vasp')
        cell = atoms.get_cell()
        symb = atoms.get_atomic_numbers()
        
        cutoffs = np.ones((len(atoms))) * 1.35  # radius around each atom (half the max bondlength)
        nl = NeighborList(cutoffs, 0., self_interaction = False, bothways = False)
        nl.update(atoms)
        
        std_len = 2.35805097505
        
        # Loop over atoms
        bond_no = 0
        bonds = np.zeros((nl.nneighbors,4))
        
        for atom1idx, (indices, offsets, atom1pos) in enumerate(zip(nl.neighbors, nl.displacements, nl.positions)):
            # Loop over bonds
            for atom2idx, offset in zip(indices, offsets):
                #if symb[atom1idx] == symb[atom2idx] == 'Si': 
                atom2pos = nl.positions[atom2idx] + np.dot(offset, cell)
                #print "pos[%i] = %.2f %.2f %.2f" % (atom1idx, atom1pos[0],atom1pos[1],atom1pos[2])
                #print "pos[%i] = %.2f %.2f %.2f" % (atom2idx, atom2pos[0],atom2pos[1],atom2pos[2])
                bond = atom2pos - atom1pos
                bond_center = atom1pos + bond/2.
                bond_z = bond_center[2]
                bondlength = np.sqrt(np.dot(bond,bond))    
                bonds[bond_no] = [bond_z, bondlength, symb[atom1idx], symb[atom2idx]]
                bond_no += 1
        self.bonds = bonds
        print bonds
    
    def get_bonds(self, atom1, atom2):
        i = np.any(np.column_stack((np.all(self.bonds[:,2:]==[atom1,atom2],axis=1), np.all(self.bonds[:,2:]==[atom2,atom1],axis=1))),axis=1)
        return self.bonds[i,0:2] # don't return the last two columns (atom numbers) since they are given


import matplotlib.pyplot as plt

bl = BondList('CONTCAR')

SiSi = bl.get_bonds(14,14)
SiSi_mean = np.mean(SiSi[:,1])
plt.axhline(SiSi_mean, color='0.60', zorder=0)
plt.axhline(2.37, color='0.70', zorder=-1, linestyle='--')
plt.scatter(SiSi[:,0], SiSi[:,1], color='0.30', zorder=2)

SiN = bl.get_bonds(7,14)
SiN_mean = np.mean(SiN[:,1])
plt.axhline(SiN_mean, color='0.60', zorder=1)
plt.axhline(1.73, color='0.70', zorder=-2, linestyle='--')
plt.scatter(SiN[:,0], SiN[:,1], color='blue', zorder=3, marker='^')

print 'Si-Si mean bond length: %.3f ' % (SiSi_mean)
print 'Si-N mean bond length: %.3f ' % (SiN_mean)


plt.ylabel(u'Bond length')
plt.xlabel(u'z [Å]')
plt.ylim(0.9*np.min(SiN[:,1]), 1.1*np.max(SiSi[:,1]) )
plt.show()
