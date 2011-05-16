#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import os,copy,sys
from copy import copy

from oppvasp.plotutils import prepare_canvas, symmetric_running_mean, symmetric_running_median
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scitools.std import seq
from oppvasp.md import Trajectory
from oppvasp.vasp.parsers import VasprunParser, IterativeVasprunParser, PoscarParser
import time,datetime
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
     RotatingMarker, ReverseBar, SimpleProgress

import psutil 

def print_memory_usage():
    p = psutil.Process(os.getpid())
    rss,vms = p.get_memory_info()
    print "Memory usage is %.1f MB" % (rss/1024.**2)

def norm(v):
    return np.sqrt(np.dot(v,v))


#############################################################################
# (1) Extract data

wdir = './'

basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
fig_filename = wdir + basename + '.pdf'
xml = wdir + 'vasprun.xml'
npz = wdir + 'trajectory_pbc.npz'

if os.path.isfile(npz):
    traj = Trajectory(filename = npz)
else:
    p = IterativeVasprunParser(xml)
    traj = p.get_all_trajectories()
    traj.save(npz)

fac = 8 * 2 * np.pi
nsteps = traj.length
natoms = traj.num_atoms

pbar = ProgressBar(widgets=['Et ögonblick...',Percentage()], maxval = nsteps).start()

lambda_x = np.zeros((nsteps))
lambda_y = np.zeros((nsteps))
lambda_z = np.zeros((nsteps))

# for each MD step
for stepno in range(nsteps):
    pbar.update(stepno)
    
    x = traj.positions[stepno,:,0]   # x coordinate of all atoms
    y = traj.positions[stepno,:,1]   # y coordinate of all atoms
    z = traj.positions[stepno,:,2]   # z coordinate of all atoms
        
    lambda_x[stepno] += np.sum(np.cos(fac * x))
    lambda_y[stepno] += np.sum(np.cos(fac * y))
    lambda_z[stepno] += np.sum(np.cos(fac * z))


pbar.finish()


lambda_x /= natoms
lambda_y /= natoms
lambda_z /= natoms
lambda_t =  1./3 * (lambda_x + lambda_y + lambda_z) 
print "Initial lattice integrity:",lambda_t[0]

#############################################################################
# (2) Plot

styles = [
       { 'color': 'black' }
   ]

prepare_canvas('10 cm') 
fig = plt.figure()

# Main plot:
p = [0.15, 0.15, 0.05, 0.05]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
ax1 = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
ax1.grid(True, which='major', color = 'gray', alpha = 0.5)
ax1.set_ylim((-.2,1.0))
ax1.set_xlabel('Time [ps]')
ax1.set_ylabel('$\lambda$')
ax1.plot(traj.time[0:nsteps]/1.e3, lambda_t, **styles[0])

#ax1.legend(('300 K', '800 K', '1300 K', '1800 K', '2300 K'), ncol = 1, loc = 'lower left', frameon = False)

plt.show()
#############################################################################
## (3) Save PDF 

sys.stdout.write("Writing %s... " % os.path.basename(fig_filename))
sys.stdout.flush()
plt.savefig(fig_filename)
sys.stdout.write("done!\n")


#print
 
