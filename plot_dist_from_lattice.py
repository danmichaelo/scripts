#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import os, sys
import numpy as np
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
     RotatingMarker, ReverseBar, SimpleProgress

from oppvasp.plotutils import prepare_canvas, symmetric_running_mean, symmetric_running_median
import matplotlib.pyplot as plt
from scipy import interpolate
from scitools.std import seq

from oppvasp.md import Trajectory
from oppvasp.vasp.parsers import VasprunParser, IterativeVasprunParser, PoscarParser

import time,datetime
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
     RotatingMarker, ReverseBar, SimpleProgress



class PlotDistanceFromLattice:
    """
    This class calculates the distances between any atom in the trajectory and the n nearest lattice site(s).
    The lattice sites are read from the file 'lattice_poscar', and the number of sites defined in the poscar
    need not be equal to the number of atoms found in the trajectory file.
    """

    def __init__(self, trajectory_dir = './', lattice_poscar = '', follow_atom = 0, step_size = 50):
        self.traj_dir = trajectory_dir
        self.lattice_poscar = lattice_poscar
        self.traj = self.read_trajectory( traj_dir = trajectory_dir )

        try:
            "%d" % follow_atom
        except TypeError:
            print "Generating geometric mean"
            geo = self.traj.get_geometric_center((0,1))
            self.traj.remove_atom(0)
            self.traj.remove_atom(0)
            follow_atom = self.traj.add_atom(geo)

        self.pos = self.traj.get_all_trajectories( coordinates = 'direct' )
        self.natoms = self.pos.shape[1]
        self.nsteps = self.pos.shape[0]

        print "Generating symmetric running mean..."
        self.pos[:,follow_atom] = symmetric_running_mean(self.pos[:,follow_atom],250)
        #pos[:,0] = symmetric_running_mean(pos[:,0],500)
        

        # Read reference lattice:
        pp = PoscarParser( lattice_poscar )
        self.lattice = pp.get_positions( coordinates = 'direct')
        self.nsites = self.lattice.shape[0]
        print "Reference POSCAR contains %d lattice sites" % (self.nsites)


        #import profile
        #profile.run('get_occupancies()')
        
        # Calculate occupancies, nearest site for atom x and distance to nearest site(s)
        self.step_size = step_size
        self.occupancies, self.followed_atom_site, self.followed_atom_r2 = self.get_occupancies( follow_atom = follow_atom)
        

    def read_trajectory(self, traj_dir = './', xml_file ='vasprun.xml', npz_pbc_file = 'trajectory_pbc.npz', npz_file = 'trajectory_nopbc.npz' ):
        if os.path.isfile(traj_dir + npz_pbc_file):
            traj = Trajectory(filename = traj_dir +npz_pbc_file)
        else:
            p = IterativeVasprunParser(traj_dir + xml_file)
            traj = p.get_all_trajectories()
            traj.save(traj_dir + npz_pbc_file)
            # we do NOT unwrap the PBCs
        return traj

    def get_occupancies(self, follow_atom = 0, closest_sites = 8):
        """
        This method generates
         (1) A (time, nsites) array with the occupancy of each lattice site at each time. 
             The occupancy of a given site is defined as the number of atoms whose positions are 
             closer to this site than any other (periodic boundary conditions are taken into account).

         (2) A (time,closest_sites) array of a given number of closest sites for atom 'follow_atom'
             at any time, and the distances to those sites.

        Parameters:
            follow_atom : (int) the atom id to follow
            closest_sites : (int) the number of closest sites to include. 
            Increasing this value does not affect the processing time very much.

        Returns:
            A (occupancies, nearest_sites, nearest_sites_r2) tuple.
            'nearest_sites' is an array containing the nearest_sites (in decreasing order) for the followed atom,
            and 'nearest_sites_r2' the corresponding distances to these sites.
        """

        n = self.traj.time[::self.step_size].shape[0]

        pbar = ProgressBar(widgets=['Et ögonblick...',Percentage()], maxval = n).start()
        occupancies = np.zeros((n,self.nsites), dtype=int)

        # initialize arrays:
        nclosest = closest_sites # plot <> closest sites
        p_site = np.zeros((n,nclosest), dtype=int)
        p_dist_r2 = np.zeros((n,nclosest))
        
        # for each MD step
        for i in range(n):
            pbar.update(i)
            stepno = i * self.step_size
            
            # for each atom
            for atno in range(self.natoms):
                atpos = self.pos[stepno,atno]
                
                # Make a (nsites,3) matrix with vectors pointing from the atom to each lattice site
                dist = self.lattice - atpos
                
                # Make a (3,nsites,3) tensor with one periodic image in each of the directions +x,-x,+y,-y,+z,-z:
                dist_i = np.abs(np.array((dist-1,dist,dist+1))) 
                # and find the shortest distance:
                dist = np.min(dist_i, axis=0)

                # Make a (nsites) vector with the distances squared (r^2)
                dist_r2 = np.sum(dist**2, axis=1)
                
                # Find closest lattice site(s):
                #closest_site = np.argmin(dist_r2)
                closest_sites = np.argsort(dist_r2)

                if atno == follow_atom:
                    p_site[i] = closest_sites[0:nclosest]
                    p_dist_r2[i] = dist_r2[p_site[i]]

                # Update occupancies array:
                occupancies[i,closest_sites[0]] += 1
           
        pbar.finish()

        p_sites = { }
        for site in p_site[:,0]: # closest site
            if site in p_sites:
                p_sites[site] += 1
            else:
                p_sites[site] = 1
        invtot = 100./p_site.shape[0]
        print "-----------"
        keys = p_sites.keys()
        for k in np.argsort(p_sites.values())[::-1]:  # reverse order 
            site = keys[k]
            print "At site %d for %.1f%% of the time" % (site,p_sites[site]*invtot)
        print "-----------"



        return occupancies, p_site, p_dist_r2


#############################################################################
# (4) Analayze 



#for i in range(n):
#    s = np.argsort(occupancies[i])
#    if occupancies[i,s[0]] == 0:
#        print "Vacancy found at site",s[0],"at step",i
#        #sys.stdout.write(str(s[0])+" ")
#    if occupancies[i,s[-2]] == 2:
#        print "Double occupation at site",s[0],"at step",i



#############################################################################
# (5) Plot

    def plot_to_file(self, fig_filename = ''):
        """
        Automatic plot method. This method is not very robust, and should in general be modified 
        to highlight the points of interest.
        """
        styles = [
               { 'color': 'black' },
               { 'color': 'gray', 'alpha': .5 }, # dashes: (ink on, ink off)
               { 'color': 'gray', 'linestyle': '--', 'dashes': (6,2), 'alpha': .5 }, # dashes: (ink on, ink off)
               { 'color': 'blue' }
           ]
        prepare_canvas('10 cm')
        fig = plt.figure()
        p = [0.15, 0.15, 0.05, 0.05]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
        ax1 = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
        ax1.set_ylabel(u'Distance [Å]')
        ax1.set_xlabel('Time [ps]')

        time = self.traj.time[::self.step_size]/1.e3
        y = np.sqrt(self.followed_atom_r2[:,0:3])*10.94

        ax1.plot(time, y, **styles[0])
        for i in range(y.shape[1]):
            ax1.axhline(y[0,i], color = 'black', linestyle='dashed', alpha = 0.5)

        y_max = 3.
        y_sub = .3

        prev_site = -1
        cur_site = -1
        cur_count = 0
        n = self.occupancies.shape[0]
        for i in range(n):
            cur_site = self.followed_atom_site[i,0]
            if cur_site != prev_site:
                if cur_count*self.step_size > 1000:
                    ax1.fill([time[i-cur_count], time[i-cur_count], time[i], time[i]], [y_max,y_max-y_sub,y_max-y_sub,y_max], color = 'green', linewidth=0., alpha = 0.2)
                    ax1.text(time[i-cur_count]+(time[i]-time[i-cur_count])/2.,y_max-(y_sub/2.),'%d' % (prev_site), horizontalalignment='center', verticalalignment='center', fontsize = 8)
                ax1.axvline(time[i], color = 'red', alpha = 0.5)
                cur_count = 0
                prev_site = cur_site
            cur_count += 1
#if label_added:
#    ax1.fill([time[i-cur_count+1], time[i-cur_count], time[i], time[i]], [y_max,y_max-y_sub,y_max-y_sub,y_max], color = 'green', linewidth=0., alpha = 0.2)

        cur_count-=1
        if cur_count*self.step_size > 1000:
            ax1.fill([time[i-cur_count], time[i-cur_count], time[i], time[i]], [y_max,y_max-y_sub,y_max-y_sub,y_max], color = 'green', linewidth=0., alpha = 0.2)
            ax1.text(time[i-cur_count]+(time[i]-time[i-cur_count])/2.,y_max-(y_sub/2.),'%d' % (prev_site), horizontalalignment='center', verticalalignment='center', fontsize = 8)


        if fig_filename == '':
            fig_filename = os.path.splitext(sys.argv[0])[0] + '.pdf'
        sys.stdout.write("Writing %s... " % os.path.basename(fig_filename))
        sys.stdout.flush()
        plt.savefig(fig_filename)
        sys.stdout.write("done!\n")

        print


if __name__ == "__main__":
    poscar = '/Users/danmichael/Documents/Studier/Master/notur/hipersol/si/si64/template/POSCAR'
    dp = PlotDistanceFromLattice('./', lattice_poscar = poscar, follow_atom = 0)
    filename = os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.pdf'
    dp.plot_to_file(filename)


