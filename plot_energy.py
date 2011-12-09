#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import os,copy,sys
from copy import copy

from oppvasp.plotutils import prepare_canvas, symmetric_running_mean, symmetric_running_median
prepare_canvas('10 cm') 
import matplotlib.pyplot as plt

import numpy as np
from scipy import interpolate
from scitools.std import seq
from oppvasp.vasp.parsers import read_trajectory
#from oppvasp.vasp.parsers import VasprunParser, IterativeVasprunParser, PoscarParser
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
     RotatingMarker, ReverseBar, SimpleProgress


import psutil 

def norm(v):
    return np.sqrt(np.dot(v,v))

wdir = './'

basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
fig_filename = wdir + basename + '.pdf'
traj = read_trajectory(wdir, unwrap_pbcs = True)
poten = (traj.total_energy - traj.kinetic_energy)/traj.num_atoms
t = traj.time/1.e3  # -> ps
poten_m = symmetric_running_mean(poten, 500)

styles = [
   { 'color': 'black', 'alpha': 0.2 },
   { 'color': 'black', 'alpha': 1. }
]

fig = plt.figure()

# Main plot:
p = [0.15, 0.15, 0.05, 0.05]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
ax1 = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
ax1.grid(True, which='major', color = 'gray', alpha = 0.5)
#ax1.set_ylim((-.2,1.0))
ax1.set_xlabel('Time [ps]')
ax1.set_ylabel('Potential energy per atom [eV]')
ax1.plot(t, poten, **styles[0])
ax1.plot(t, poten_m, **styles[1])

#plt.show()

sys.stdout.write("Writing %s... " % os.path.basename(fig_filename))
sys.stdout.flush()
plt.savefig(fig_filename)
sys.stdout.write("done!\n")

