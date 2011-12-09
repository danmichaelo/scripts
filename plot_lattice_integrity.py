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


lambda_x = np.zeros((nsteps))
lambda_y = np.zeros((nsteps))
lambda_z = np.zeros((nsteps))



### <!-- Test begin --> 
from oppvasp.vasp.parsers import PoscarParser
ps = PoscarParser('/Users/danmichael/master/notur/hipersol/templates/si64/POSCAR').get_structure()
r = ps.get_positions('direct')

# Shift all positions
r += -0.01+np.random.rand(64,3)*0.02
natoms = r.shape[0]
test1 = 1./3 * (np.sum(np.cos(fac * r[:,0]))/natoms + np.sum(np.cos(fac * r[:,1]))/natoms + np.sum(np.cos(fac * r[:,2]))/natoms )

kvec = np.array((8,0,0)) # 2**3 for 2x2x2 supercell
s = 2*np.pi*np.dot(r, kvec)

test2 = np.sum(np.exp(np.complex(0,1) * s)).real / natoms

cossum = np.sum(np.cos(s))/natoms
sinsum = np.sum(np.sin(s))/natoms
test3 = np.sqrt(cossum**2+sinsum**2)

print test1
print test2
print test3



pos = traj.positions

T1 = np.zeros(pos.shape[0])
T2 = np.zeros(pos.shape[0])
T3 = np.zeros(pos.shape[0])
#r = traj.positions[-1]

pbar = ProgressBar(widgets=['Et ögonblick...',Percentage()], maxval = nsteps).start()
for idx, r in enumerate(pos):
    pbar.update(idx)
    natoms = r.shape[0]
    T1[idx] = 1./3 * (np.sum(np.cos(fac * r[:,0]))/natoms + np.sum(np.cos(fac * r[:,1]))/natoms + np.sum(np.cos(fac * r[:,2]))/natoms )

    kvec = np.array((8,0,0)) # 2**3 for 2x2x2 supercell
    s = 2*np.pi*np.dot(r, kvec)

    #print "T2:",
    T2[idx] = np.sum(np.exp(np.complex(0,1) * s)).real / natoms

    cossum = np.sum(np.cos(s))/natoms
    sinsum = np.sum(np.sin(s))/natoms
    T3[idx] = np.sqrt(cossum**2+sinsum**2)

T1m = symmetric_running_mean(T1, 500)
T2m = symmetric_running_mean(T2, 500)
T3m = symmetric_running_mean(T3, 500)

## for each MD step
#for stepno in range(nsteps):
    
#    x = traj.positions[stepno,:,0]   # x coordinate of all atoms
#    y = traj.positions[stepno,:,1]   # y coordinate of all atoms
#    z = traj.positions[stepno,:,2]   # z coordinate of all atoms
        
#    lambda_x[stepno] += np.sum(np.cos(fac * x))
#    lambda_y[stepno] += np.sum(np.cos(fac * y))
#    lambda_z[stepno] += np.sum(np.cos(fac * z))


pbar.finish()


#lambda_x /= natoms
#lambda_y /= natoms
#lambda_z /= natoms
#lambda_t =  1./3 * (lambda_x + lambda_y + lambda_z) 
#print "Initial lattice integrity:",lambda_t[0]

#############################################################################
# (2) Plot

styles = [
        { 'color': 'black', 'alpha': 0.2 },
        { 'color': 'green', 'alpha': 0.2 },
        { 'color': 'blue', 'alpha': 0.2 },
        { 'color': 'black' },
       { 'color': 'green' },
       { 'color': 'blue' }
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
ax1.plot(traj.time[0:nsteps]/1.e3, T1, **styles[0])
ax1.plot(traj.time[0:nsteps]/1.e3, T2, **styles[1])
ax1.plot(traj.time[0:nsteps]/1.e3, T3, **styles[2])
ax1.plot(traj.time[0:nsteps]/1.e3, T1m, label = 'T1', **styles[3])
ax1.plot(traj.time[0:nsteps]/1.e3, T2m, label = 'T2', **styles[4])
ax1.plot(traj.time[0:nsteps]/1.e3, T3m, label = 'T3', **styles[5])

ax1.legend( ncol = 1, frameon = False )

plt.show()
#############################################################################
## (3) Save PDF 

sys.stdout.write("Writing %s... " % os.path.basename(fig_filename))
sys.stdout.flush()
plt.savefig(fig_filename)
sys.stdout.write("done!\n")


#print
 
