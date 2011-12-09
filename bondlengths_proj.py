#!/usr/bin/env python
#encoding=utf8

import numpy as np
from oppvasp.util import get_pairs

def lensq(arr):
    """Returns the squared length of the vectors along the rows in
    *arr*, where *arr* is a 2D array."""
    return (arr**2).sum(axis=1)

def main():
    import ase
    #from connectivities import connectivities, partition, plot
    atoms = ase.io.read('CONTCAR', format='vasp')
    #dists, idx, bonds = connectivities(atoms, maxdist=2.5)
    #groups = partition(atoms, dists, idx, mincount=50, maxwidth=2.0)
    #cats = categorize(atoms, dists, idx)
    #cats[0]

    pos = atoms.positions

    # Get all atom pairs:
    n = len(atoms)
    pairs = get_pairs(n)

    # Filter pairs based on atom types:
    bond = (14,14)
    pnum = atoms.numbers[pairs]
    # numbers == (bond[0],bond[1]) or numbers == (bond[1],bond[0]) :
    mask = np.any(np.column_stack((np.all(pnum == bond, axis=1), np.all(pnum == reversed(bond), axis = 1))),axis=1)
    del pnum

    # Calculate distancess:
    neighbours = [(-1, 2), (-1, 2), (-1, 2)]
    vectors = atoms.positions[pairs[:,0]] - atoms.positions[pairs[:,1]]
    tmp = np.zeros((len(vectors), 2))
    tmp[:,0] = lensq(vectors)
    for d0 in range(*neighbours[0]):
        for d1 in range(*neighbours[1]):
            for d2 in range(*neighbours[2]):
                if d0 == d1 == d2 == 0:
                    continue
                tmp[:,1] = lensq(vectors + np.dot(np.array([d0, d1, d2]), 
                                               atoms.cell))
                tmp.min(axis=1, out=tmp[:,0])

    #
    # vectors må korrigeres for PBC. Sjekk om vi kan få pair list rett ut fra ASE
    #

    dists = np.sqrt(tmp[:,0])
    del tmp

    # Filter pairs based on bond lengths:
    maxdist = 2.8
    if maxdist is not None:
        # add to mask
        mask = np.all(np.column_stack((mask, dists <= maxdist)), axis = 1)

    dists = dists[mask]
    pairs = pairs[mask,:]
    vectors = vectors[mask,:]

    import matplotlib.pyplot as plt

    vec2d = np.zeros((vectors.shape[0],2))
    x = np.array([1,0,0])
    y = np.array([0,1,0])
    vec2d = vectors[:,::2]

    pos2d = pos[pairs[:,0],::2]

    print vec2d.shape
    print pos2d.shape

    plt.quiver(pos[:,0], pos[:,1], vec2d[:,0], vec2d[:,1])

    plt.show()


if __name__ == '__main__':
    main()
    
