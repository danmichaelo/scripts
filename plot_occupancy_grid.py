#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import os, sys
import numpy as np
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
     RotatingMarker, ReverseBar, SimpleProgress

from oppvasp.plotutils import prepare_canvas, symmetric_running_mean, symmetric_running_median
from matplotlib import pyplot as plt, mpl
from scipy import interpolate
from scitools.std import seq

from oppvasp.md import Trajectory
from oppvasp.vasp.parsers import VasprunParser, IterativeVasprunParser, PoscarParser

import time,datetime
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
     RotatingMarker, ReverseBar, SimpleProgress

class PlotOccupancyGrid:

    def __init__(self, trajectory_dir = './', lattice_poscar = '', follow_atom = 0, step_size = 200):

        # Read reference lattice:
        pp = PoscarParser( lattice_poscar )
        self.lattice = pp.get_positions( coordinates = 'direct')
        self.nsites = self.lattice.shape[0]
        print "Reference POSCAR contains %d lattice sites" % (self.nsites)
        
        # Get occupancies
        self.traj_dir = trajectory_dir
        self.lattice_poscar = lattice_poscar
        self.traj = self.read_trajectory( traj_dir = trajectory_dir )
        self.occ, self.inhabitants = self.traj.get_occupancies( lattice = self.lattice, step_size = step_size )
        self.time = self.traj.time[::step_size]/1.e3
        #self.vel = self.traj.get_velocities()

        #self.pos = self.traj.get_all_trajectories( coordinates = 'direct' )

        #print "Generating symmetric running mean..."
        #self.pos[:,follow_atom] = symmetric_running_mean(self.pos[:,follow_atom],250)
        #pos[:,0] = symmetric_running_mean(pos[:,0],500)
        
        #import profile
        #profile.run('get_occupancies()')
        
        # Calculate occupancies, nearest site for atom x and distance to nearest site(s)
        #self.step_size = step_size
        #self.occupancies, self.followed_atom_site, self.followed_atom_r2 = self.get_occupancies( follow_atom = follow_atom)
        

    def read_trajectory(self, traj_dir = './', xml_file ='vasprun.xml', npz_pbc_file = 'trajectory_pbc.npz', npz_file = 'trajectory_nopbc.npz' ):
        if os.path.isfile(traj_dir + npz_pbc_file):
            traj = Trajectory(filename = traj_dir +npz_pbc_file)
        else:
            p = IterativeVasprunParser(traj_dir + xml_file)
            traj = p.get_all_trajectories()
            traj.save(traj_dir + npz_pbc_file)
            # we do NOT unwrap the PBCs
        return traj


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
        
        #
        #ax1.set_ylabel(u'Distance [Å]')
        #ax1.set_xlabel('Time [ps]')
        #print self.occ.shape
        #print self.occ
        #plt.pcolor(C, cmap=mpl.cm.gray_r )
        # http://matplotlib.sourceforge.net/examples/pylab_examples/pcolor_demo2.html
        #p = [0.15, 0.15, 0.05, 0.52]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
        #ax1 = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
        #norm = mpl.colors.Normalize(vmin=0, vmax=3)
        #im = ax1.imshow(self.occ.T, cmap=mpl.cm.gray_r, norm=norm, aspect='auto', interpolation='nearest')

        nsteps = self.occ.shape[0]
        nsites = self.occ.shape[1]
        natoms = self.traj.num_atoms

        lw = 10.
        
        #
        p = [0.15, 0.05, 0.05, 0.05]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
        ax2 = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
        #print self.vel.shape
        #ax2.acorr(self.vel[:,0,0],usevlines=True, lw=2)
        #norm2 = mpl.colors.Normalize(vmin=0, vmax=natoms)
        
        maxocc = 1
        imgrid = np.ones((nsteps, nsites*maxocc*lw,3)) # RGB float32 array

        colors=np.ones((natoms,3))*.9   # default atom color
        colors[0] = [1.,0.,0.]          # color for atom 0


        print self.occ.shape
        print self.inhabitants.shape
        print imgrid.shape
        for i in range(self.occ.shape[0]):
            for j in range(self.occ.shape[1]):
                #imgrid[i,j*maxocc*lw:j*maxocc*lw+lw] = colors[self.inhabitants[i,j,0]]
                for k in range(self.occ[i,j]):
                    imgrid[i,(j*maxocc+k)*lw:(j*maxocc+k)*lw+lw] = colors[self.inhabitants[i,j,k]] # 0,natoms -> 1,natoms+1 so 0=no atoms
                    if k+1 >= maxocc:
                        print "Warning: Too high occupancy"
                        break
                #    #print (j*4+k)/4., colors[self.inhabitants[i,j,k]]
        im2 = ax2.imshow(np.transpose(imgrid, axes=[1,0,2]), aspect='auto', interpolation='nearest')
        #im2.
        #im2 = ax2.imshow(np.transpose(imgrid, axes=[1,0,2]), aspect='auto', interpolation='nearest')
        #ax2.set_ylim(-1*maxocc*lw,nsites*maxocc*lw+maxocc+lw)
        #ax2.minorticks_on()
        #ax2.set_yticks(np.arange(0,64*maxocc*lw,4), minor=True)
        #ax2.set_yticks(np.arange(0,64*4,10*4))
        #ax2.set_yticks([])
        #ax2.minorticks_off()
        #ax1.


        #ticks = ax2.get_yticks()
        #print ticks
        #labels = [ t/4. for t in ticks ]
        #print labels
        #ax2.set_yticklabels(labels)

        

#        y_max = 3.
#        y_sub = .3

#        prev_site = -1
#        cur_site = -1
#        cur_count = 0
#        n = self.occupancies.shape[0]
#        for i in range(n):
#            cur_site = self.followed_atom_site[i,0]
#            if cur_site != prev_site:
#                if cur_count*self.step_size > 1000:
#                    ax1.fill([time[i-cur_count], time[i-cur_count], time[i], time[i]], [y_max,y_max-y_sub,y_max-y_sub,y_max], color = 'green', linewidth=0., alpha = 0.2)
#                    ax1.text(time[i-cur_count]+(time[i]-time[i-cur_count])/2.,y_max-(y_sub/2.),'%d' % (prev_site), horizontalalignment='center', verticalalignment='center', fontsize = 8)
#                ax1.axvline(time[i], color = 'red', alpha = 0.5)
#                cur_count = 0
#                prev_site = cur_site
#            cur_count += 1
##if label_added:
##    ax1.fill([time[i-cur_count+1], time[i-cur_count], time[i], time[i]], [y_max,y_max-y_sub,y_max-y_sub,y_max], color = 'green', linewidth=0., alpha = 0.2)

#        cur_count-=1
#        if cur_count*self.step_size > 1000:
#            ax1.fill([time[i-cur_count], time[i-cur_count], time[i], time[i]], [y_max,y_max-y_sub,y_max-y_sub,y_max], color = 'green', linewidth=0., alpha = 0.2)
#            ax1.text(time[i-cur_count]+(time[i]-time[i-cur_count])/2.,y_max-(y_sub/2.),'%d' % (prev_site), horizontalalignment='center', verticalalignment='center', fontsize = 8)


        if fig_filename == '':
            fig_filename = os.path.splitext(sys.argv[0])[0] + '.pdf'
        sys.stdout.write("Writing %s... " % os.path.basename(fig_filename))
        sys.stdout.flush()
        plt.savefig(fig_filename)
        sys.stdout.write("done!\n")

        print


if __name__ == "__main__":
    poscar = '/Users/danmichael/Documents/Studier/Master/notur/hipersol/si/si64/template/POSCAR'
    dp = PlotOccupancyGrid('./', lattice_poscar = poscar)
    filename = os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.pdf'
    dp.plot_to_file(filename)


