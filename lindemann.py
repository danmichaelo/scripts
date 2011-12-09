#!/usr/bin/env python
#encoding=utf8

import numpy as np
from oppvasp.vasp.parsers import read_trajectory
from oppvasp.util import get_pairs

from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
     RotatingMarker, ReverseBar, SimpleProgress

def lensq(arr):
    """Returns the squared length of the vectors along the rows in
    *arr*, where *arr* is a 2D array."""
    return (arr**2).sum(axis=1)

def lindemann_index(trajectory_dir='./'):
    """
    Calculates the Lindemann index
    \delta = \frac{2}{N(N-1)} \sum_{i<j} \frac{\sqrt{\langle r_{ij}^2\rangle_t - \langle r_{ij} \rangle^2_t }}{ \langle r_{ij} \rangle_t }
    """
    
    # Import vasprun.xml:
    traj = read_trajectory( dir = trajectory_dir, unwrap_pbcs = False )
    pos = traj.get_all_trajectories( coords = 'direct' )
    
    nsteps = pos.shape[0]
    natoms = pos.shape[1]

    # Get all atom pairs:
    pairs = get_pairs(natoms)
    npairs = pairs.shape[0]

    print "%d steps, %d atoms, %d pairs" % (nsteps, natoms, npairs)

    stepsize = 50
    samplesteps = np.arange(0,nsteps,stepsize)
    r = np.zeros((samplesteps.size, npairs)) # increase stepsize if we run out of memory
    r2 = np.zeros((samplesteps.size, npairs)) # increase stepsize if we run out of memory
    print "Allocated %.3f MB" % ((r.nbytes+r2.nbytes)/(1024.**2))

    pbar = ProgressBar(widgets=['At step...',SimpleProgress()], maxval = samplesteps.size).start()
    for idx,s in enumerate(samplesteps):
        pbar.update(idx)
        x = pos[s,pairs[:,0]] - pos[s,pairs[:,1]]
        x = x - (2*x).astype('int')  # minimum image convention
        r2[idx] = lensq(x)
        r[idx] = np.sqrt(r2[idx])
    pbar.finish()

    meanr2 = np.mean(r2,axis=0)
    meanr = np.mean(r,axis=0)

    return np.sum(np.sqrt((meanr2 - meanr**2))/meanr)*2/(npairs*2)


    #print np.mean(r2)
    # Filter pairs based on atom types:
    #bond = (14,14)
    #pnum = atoms.numbers[pairs]
    ## numbers == (bond[0],bond[1]) or numbers == (bond[1],bond[0]) :
    #mask = np.any(np.column_stack((np.all(pnum == bond, axis=1), np.all(pnum == reversed(bond), axis = 1))),axis=1)
    #del pnum

    # Calculate \vec{r}_{ij} and r_{ij}^2



    #dists = np.sqrt(tmp[:,0])

    # Filter pairs based on bond lengths:
    #maxdist = 2.8
    #if maxdist is not None:
    #    # add to mask
    #    mask = np.all(np.column_stack((mask, dists <= maxdist)), axis = 1)

    #dists = dists[mask]
    #pairs = pairs[mask,:]

    #print np.mean(dists)
    #import matplotlib.pyplot as plt



if __name__ == '__main__':
    print lindemann_index('./')
    
