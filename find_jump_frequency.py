#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import os, sys
import numpy as np
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
     RotatingMarker, ReverseBar, SimpleProgress

import matplotlib.pyplot as plt
#from oppvasp.plotutils import prepare_canvas, symmetric_running_mean, symmetric_running_median
from oppvasp.vasp.parsers import read_trajectory

from scipy import interpolate
from scitools.std import seq

from oppvasp.md import Trajectory
from oppvasp.vasp.parsers import VasprunParser, IterativeVasprunParser, PoscarParser

import time,datetime

from scipy import fft

class FindNearestSite(object):
    """
    This class calculates the distances between any atom in the trajectory and the n nearest lattice site(s).
    The lattice sites are read from the file 'lattice_poscar', and the number of sites defined in the poscar
    need not be equal to the number of atoms found in the trajectory file.
    """

    def __init__(self, trajectory_dir = './', lattice_poscar = '', atomid = 0):
        """
        time_steps in ps
        """
        self.traj_dir = trajectory_dir
        self.lattice_poscar = lattice_poscar
        self.traj = read_trajectory( dir = trajectory_dir, unwrap_pbcs = True )

        self.time = self.traj.time/1000.
        self.pos = self.traj.get_all_trajectories( coords = 'cart' )
        self.natoms = self.pos.shape[1]
        self.nsteps = self.pos.shape[0]

        #print "Generating symmetric running mean..."
        #self.pos[:,atomid] = symmetric_running_mean(self.pos[:,atomid],250)
        #pos[:,0] = symmetric_running_mean(pos[:,0],500)
        

        # Read reference lattice:
        pp = PoscarParser( lattice_poscar )
        self.lattice = pp.get_positions( coords = 'cart')
        self.nsites = self.lattice.shape[0]
        #print "Reference POSCAR contains %d lattice sites" % (self.nsites)

        #import profile
        #profile.run('get_occupancies()')
        
        # Calculate occupancies, nearest site for atom x and distance to nearest site(s)
        self.nearest_site, self.nearest_site_r = self.find_nearest_sites( atomid = atomid)

    def find_vib_freq(self, atomid = 0, firstframe = 1000, lastframe = 10000):

        time = self.time[firstframe:lastframe+1]
        #(np.arange(lastframe-firstframe+1)+firstframe)*self.timestep
        atompos = self.pos[firstframe:lastframe+1,atomid]
        siteid, sitedist = self.find_nearest_site(atomid, 0)
        print "Nearest site is %d (%.3f Ang)" % (siteid, sitedist)
        sitepos = self.lattice[siteid]
        dist = np.sqrt(np.sum((sitepos-atompos)**2, axis=1))
        print dist.shape
        print "Plotting vib"

        plt.subplot(2,1,1)
        plt.plot(time,dist)
        plt.xlim(time[0],time[-1])
        #print np.linspace(time[0], time[-1], 21)
        #plt.xticks(np.linspace(time[0], time[-1], 10))
        plt.xlabel("Time [ps]")
        plt.ylabel(u"$r(t)-R$ [Å]")
        plt.grid('on')
        #xlabels = plt.get_xticklabels() 
        #for label in xlabels: 
        #    label.set_rotation(45) 
        
        #n = dist.size
        #w = np.fft.fft(dist)
        #w = np.fft.fftshift(w)
        #freq = np.fft.fftfreq(n, d=self.timestep)

        qfft = fft(dist)
        n=len(qfft)
        print "Len is",n
        power = np.abs(qfft)**2
        
        #nyquist=1.0/2
        #freq=np.array(range(n/2))/(n/2.0)*nyquist
        freq = np.fft.fftfreq(n)

        #power = power[1:]
        #freq = freq[1:]

        #period=1./freq
        
        #maxidx = np.argmax(power)
        #print "Strongest period:",period[maxidx]
        

        #for coef,freq in zip(w,freqs):
        #    if coef:
        #        print('{c:>6} * exp(2 pi i t * {f})'.format(c=coef,f=freq))

        #df = rfft(dist)
        plt.subplot(2,1,2)
        #plt.hist(power[1:n/2], 50)
        plt.plot(freq[1:n/2], power[1:n/2], '-b')
        #plt.plot(freq, qfft.imag, '-g')
        plt.xlim(0.0, 0.01)
        plt.savefig('out.pdf')
        return dist
    
    def find_nearest_site(self, atomid = 0, frameid = 0):
        atpos = self.pos[frameid,atomid]
        
        # Make a (nsites,3) matrix with vectors pointing from the atom to each lattice site
        dist = self.lattice - atpos

        # Make a (nsites) vector with the distances squared (r^2)
        dist_r = np.sqrt(np.sum(dist**2, axis=1))
        
        # Find closest lattice site(s):
        closest_site = np.argmin(dist_r)

        return closest_site, dist_r[closest_site]

    def find_nearest_sites(self, atomid = 0):
        n = self.traj.time.shape[0]
        closest_site = np.zeros((n), dtype=int)
        closest_site_r = np.zeros((n))
        for i in range(n):
            closest_site[i], closest_site_r[i] = self.find_nearest_site(atomid, i)
           
        return closest_site, closest_site_r

    
    def analyze(self, min_time, max_dist = 1., verbose = True):
        """
        This function will count the number of "jumps" an atom makes between lattice sites.
        The lattice sites themselves are taken to be completely static, and possible drift 
        of the whole lattice will not be monitored here, but should be checked by other means.

        Definitions
        ----------
        * An atom A _sits_ at lattice site X when X is the nearest lattic site to A
          for longer than <min_time>, and (optionally) the distance between A and X 
          has been less than <max_dist> during that time.

        * An atom _jumps_ when it changes from _sitting_ at one site to _sitting_ at 
          another site.

        Parameters
        ----------
        min_time : float
            the minimum time the atom has to stay near a site to be defined as sitting on that site (in picoseconds)
        max_dist : float
            the maxmimum distance an atom can have from a site to be defined as sitting on that site (in Angstrom)
            Can be set to none to avoid this criterion.
        verbose : bool
            print information about every single jump if True
        """
        nsteps = self.nearest_site.shape[0]
        cur_site = -1
        prev_site = -1
        stepcount = 0
        jumpcount = 0
        for i in range(nsteps):
            cur_site = self.nearest_site[i]
            if cur_site != prev_site and prev_site != -1:
                # looks like a jump
                t = self.time[i] - self.time[i-stepcount] # could be made simpler, but allows for variable timestep sizes (which is probably never used though)
                rmin = np.min(self.nearest_site_r[i-stepcount:i])
                # is it really a jump?
                if t > min_time and (max_dist == None or rmin < max_dist):
                    if verbose:
                        print u"Jumping to site %d after a stay of %.3f ps (%d steps) at site %s (closest encounter: %.2f Å)" %(cur_site, t, stepcount, prev_site, rmin)
                    jumpcount += 1
                stepcount = 0
            prev_site = cur_site
            stepcount += 1
        tottime = self.time[-1] + 0.001  # since we start at zero
        return jumpcount, tottime
        
if __name__ == "__main__":
    poscar = '/Users/danmichael/Documents/Studier/Master/notur/hipersol/templates/si64/POSCAR_3x3_symmetric'
    dp = FindNearestSite('./', lattice_poscar = poscar, atomid = 0)
    d = dp.find_vib_freq(32, firstframe = 0, lastframe = 29999)
    plt.show()
    #dp.analyze( min_time = 0.5)

    #filename = os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.pdf'
    #dp.plot_to_file(filename)


