#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

import os, sys
import numpy as np
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
     RotatingMarker, ReverseBar, SimpleProgress

from oppvasp.plotutils import prepare_canvas, symmetric_running_mean, symmetric_running_median
from oppvasp.vasp.parsers import read_trajectory

from scipy import interpolate
from scitools.std import seq

from oppvasp.md import Trajectory
from oppvasp.vasp.parsers import VasprunParser, IterativeVasprunParser, PoscarParser

import time,datetime

class DistanceFromLatticeSites(object):
    """
    This class calculates the distances between any atom in the trajectory and the n nearest lattice site(s).
    The lattice sites are read from the file 'lattice_poscar', and the number of sites defined in the poscar
    need not be equal to the number of atoms found in the trajectory file.
    """

    def __init__(self, trajectory_dir = './', lattice_poscar = ''):
        
        
        self.traj_dir = trajectory_dir
        self.lattice_poscar = lattice_poscar
        self.traj = read_trajectory( trajectory_dir, unwrap_pbcs = True )

        self.pos = self.traj.get_all_trajectories( coords = 'cartesian' )
        self.natoms = self.pos.shape[1]
        self.nsteps = self.pos.shape[0]

        self.poscar = PoscarParser( lattice_poscar )
        
        # Read reference lattice:
        self.lattice = self.poscar.get_positions( coords = 'cartesian')
        self.nsites = self.lattice.shape[0]
        print "Reference POSCAR contains %d lattice sites" % (self.nsites)

        self.bottombar_height = .3
        self.bottombar_colors = {
            892: 'green',
            860: 'yellow',
            868: 'blue',
            844: 'red'
        }
    
    def read_trajectory(self, traj_dir = './', xml_file ='vasprun.xml', npz_pbc_file = 'trajectory_pbc.npz', npz_file = 'trajectory_nopbc.npz' ):
        if os.path.isfile(traj_dir + npz_nopbc_file):
            traj = Trajectory(filename = traj_dir +npz_nopbc_file)
        else:
            p = IterativeVasprunParser(traj_dir + xml_file)
            traj = p.get_all_trajectories( coords = 'cart' )
            traj.save(traj_dir + npz_pbc_file)
            # we do NOT unwrap the PBCs
        return traj

    def follow(self, follow_atom = 0, num_neighbours = 3, step_size = 50):
        self.num_neighbours = num_neighbours # numbers of nearest neighbours to plot distance to

        #try:
        #    "%d" % follow_atom
        #except TypeError:
        #    print "Generating geometric mean"
        #    geo = self.traj.get_geometric_center((0,1))
        #    self.traj.remove_atom(0)
        #    self.traj.remove_atom(0)
        #    follow_atom = self.traj.add_atom(geo)

        print "Generating symmetric running mean..."
        self.pos[:,follow_atom] = symmetric_running_mean(self.pos[:,follow_atom],250)
        #pos[:,0] = symmetric_running_mean(pos[:,0],500)
        
        #import profile
        #profile.run('get_occupancies()')
        
        # Calculate occupancies, nearest site for atom x and distance to nearest site(s)
        self.step_size = step_size
        self.occupancies, self.followed_atom_site, self.followed_atom_r2 = self.get_occupancies( follow_atom = follow_atom)

    def get_index(self,lst,val):
        """Returns first index of val in lst"""
        return [i for i,x in enumerate(lst) if x == val][0]

    def get_occupancies(self, follow_atom = 0, closest_sites = 8, make_periodic_supercell = False):
        """
        This method generates
         (1) A (time, nsites) array with the occupancy of each lattice site at each time. 
             The occupancy of a given site is defined as the number of atoms whose positions are 
             closer to this site than any other (periodic boundary conditions are taken into account).

         (2) A (time,closest_sites) array of a given number of closest sites for atom 'follow_atom'
             at any time, and the distances to those sites.

        Parameters:
            follow_atom : (int) or list of ints - the atom id(s) to follow
            closest_sites : (int) the number of closest sites to include. 
            Increasing this value does not affect the processing time very much.

        Returns:
            A (occupancies, nearest_sites, nearest_sites_r2) tuple.
            'nearest_sites' is an array containing the nearest_sites (in decreasing order) for the followed atom,
            and 'nearest_sites_r2' the corresponding distances to these sites.
        """

        n = self.traj.time[::self.step_size].shape[0]

        pbar = ProgressBar(widgets=['Frame...',SimpleProgress(),' (step size: %d)' % (self.step_size)], maxval = n*self.step_size).start()
        occupancies = np.zeros((n,self.nsites), dtype=int)

        # initialize arrays:
        nclosest = closest_sites # plot <> closest sites
        p_site = np.zeros((len(follow_atom),n,nclosest), dtype=int)
        p_dist_r2 = np.zeros((len(follow_atom),n,nclosest))
        # for each MD step
        for i in range(n):
            stepno = i * self.step_size
            pbar.update(stepno)

            
            # for each atom
            for atno in range(self.natoms):
                dist = self.lattice - self.pos[stepno,atno]
                atpos = self.pos[stepno,atno]
                
                # Make a (nsites,3) matrix with vectors pointing from the atom to each lattice site
                dist = self.lattice - atpos

                # a) Minimum image convention work for most unit cells:
                #    (use direct coordinates)
                #
                #dist = dist - (2*dist-1).astype(np.int)
               
                # b) Alternative (slower) routine that can be used with very skewed unit cells: 
                #    (use cartesian instead of direct coordinates)
                #
                if make_periodic_supercell:
                    # Make a (3,nsites,3) tensor with one periodic image in each of the directions +x,-x,+y,-y,+z,-z:
                    dist_i = np.abs(np.array((dist-1,dist,dist+1))) 
                    # and find the shortest distance:
                    dist = np.min(dist_i, axis=0)

                # Make a (nsites) vector with the distances squared (r^2)
                dist_r2 = np.sum(dist**2, axis=1)
                
                # Find closest lattice site(s):
                #closest_site = np.argmin(dist_r2)
                closest_sites = np.argsort(dist_r2)

                if atno in follow_atom:
                    idx = self.get_index(follow_atom,atno)
                    p_site[idx,i] = closest_sites[0:nclosest]
                    p_dist_r2[idx,i] = dist_r2[p_site[idx,i]]

                # Update occupancies array:
                occupancies[i,closest_sites[0]] += 1
           
        pbar.finish()
        # (atno, stepno, nclosest)

        #p_sites = { }
        ## loop over atno
        #for atidx, at in enumerate(p_site[:,:,0]):
        #    for site in at[:,0]: # closest site
        #        if site in p_sites:
        #            p_sites[site] += 1
        #        else:
        #            p_sites[site] = 1
        #    invtot = 100./p_site.shape[0]
        #    print "-----------"
        #    keys = p_sites.keys()
        #    for k in np.argsort(p_sites.values())[::-1]:  # reverse order 
        #        site = keys[k]
        #        print "At site %d for %.1f%% of the time (%.2f,%.2f,%.2f)" % (site,p_sites[site]*invtot,self.lattice[site][0],self.lattice[site][1],self.lattice[site][2])
        #    print "-----------"


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

    def bottombar_fill(self, plt, index, t0, t, labs):
        try:
            color = self.bottombar_colors[index]
        except KeyError:
            color = 'yellow'
        #print index,t0,t
        y0 = - self.bottombar_height
        plt.fill( [t0,t0,t,t], [0,y0,y0,0], color = color, linewidth = 0., alpha = 0.3 )
        
        if index not in labs and t-t0 > 0.8:
            try:
                label = str(self.bottombar_labels[index])
            except KeyError:
                label = str(index)

            print " -> Adding label",label,"(atom %d)"%(index)
            plt.text(t0 + (t-t0)/2., y0/2., label, horizontalalignment = 'center', verticalalignment = 'center', fontsize = 8)
            labs.append(index)
        return labs

    def plot(self):
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

        p = [0.15, 0.57, 0.05, 0.03]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
        upper = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
        upper.set_ylabel(u'$r(t)-r(0)$ [Å]')
        upper.set_xticklabels([])

        p = [0.15, 0.15, 0.05, 0.45]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
        lower = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
        lower.set_xlabel('Time $t$ [ps]')
        lower.set_ylabel(u'$r(t)-R$ (å)')

        time = self.traj.time[::self.step_size]/1.e3
        y = np.sqrt(self.followed_atom_r2[:,0:self.num_neighbours])

        lower.plot(time, y, **styles[0])
        ticks = [0.]
        ticklabels = [0.]
        for i in range(y.shape[1]):
            lower.axhline(y[0,i], color = 'black', linestyle='dashed', alpha = 0.5)
            print "Y:",y[0,i]
            ticks.append(y[0,i])
            ticklabels.append("%.2f" % (y[0,i]))
        lower.set_yticks(ticks)
        lower.set_yticklabels(ticklabels)

        r = np.sqrt(np.sum((self.pos[::self.step_size,0]-self.pos[0,0])**2, axis=1))

        print time.shape
        print r.shape
        upper.plot(time, r)

        # ================== Paint bottombar ==================================

        prev_site = -1
        cur_site = -1
        cur_count = 0
        n = self.occupancies.shape[0]
        labelled_sites = []
        for i in range(n):
            cur_site = self.followed_atom_site[i,0]
            if cur_site != prev_site and prev_site != -1:
                labelled_sites = self.bottombar_fill(lower, prev_site, time[i-cur_count], time[i], labelled_sites)
                cur_count = 0
            prev_site = cur_site
            cur_count += 1
        
        #cur_count-=1
        labelled_sites = self.bottombar_fill(lower, prev_site, time[-cur_count], time[-1], labelled_sites)

        self.plt = plt
        self.lower = lower
        self.upper = upper

    def save(self, fig_filename = ''):
        if fig_filename == '':
            fig_filename = os.path.splitext(sys.argv[0])[0] + '.pdf'
        sys.stdout.write("Writing %s... " % os.path.basename(fig_filename))
        sys.stdout.flush()
        self.plt.savefig(fig_filename)
        sys.stdout.write("done!\n")

    def time_at_lattice_sites(self, treshold = 0.5):
        print self.followed_atom_r2.shape  # (65, 600, 8)
        y = np.sqrt(self.followed_atom_r2[:,:,0:self.num_neighbours])
        print "y",y.shape
        return [float(j[j<treshold].shape[0])/j.shape[0] for j in y]
        #print y[not y>treshold].shape


if __name__ == "__main__":
    poscar = '/Users/danmichael/Documents/Studier/Master/notur/hipersol/templates/si64/POSCAR_3x3_symmetric'
    dp = DistanceFromLatticeSites('./', lattice_poscar = poscar, follow_atom = 0, step_size = 50, num_neighbours = 1)
    filename = os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.pdf'
    dp.plot()
    dp.save(filename)
    #dp.analyse()

