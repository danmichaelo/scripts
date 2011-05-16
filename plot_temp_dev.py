#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

from oppvasp import plotutils
import os,copy,sys
from copy import copy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scitools.std import seq
import time,datetime

import psutil 

from oppvasp import getAtomicNumberFromSymbol
from oppvasp.vasp.parsers import IterativeVasprunParser, PoscarParser
from oppvasp.md import Trajectory

def norm(v):
    return np.sqrt(np.dot(v,v))

#############################################################################
# (1) Extract data

d = './'
xml_filename = d + 'vasprun.xml'
npz_filename = d + 'trajectory_pbc.npz'
fig_filename = os.path.splitext(sys.argv[0])[0] + '.pdf'

if os.path.isfile(npz_filename):
    traj = Trajectory(filename = npz_filename)
else:
    p = IterativeVasprunParser(xml_filename)
    traj = p.get_all_trajectories()
    traj.save(npz_filename)

basis = traj.basis
atoms = traj.atoms
nsteps = traj.length # NSW
nions = traj.num_atoms
positions = traj.positions
etot = traj.total_energy
ekin = traj.kinetic_energy
time = traj.time

# min/max x:
r_min = 0
r_max = 5000


#############################################################################
# (2) Plot

styles = [
        { 'color': 'gray', 'alpha': .5 },
        { 'color': 'black', 'linestyle': '--', 'dashes': (4,1), 'alpha': .8 }, # dashes: (ink on, ink off)
        { 'color': 'black' },
        { 'color': 'gray', 'alpha': .5 } # dashes: (ink on, ink off)
    ]

plotutils.prepare_canvas( width = '7 cm', fontsize = 9, fontsize_small = 8 ) 
fig = plt.figure()

toten = traj.total_energy[r_min:r_max]
toten_ra = plotutils.symmetric_running_mean(toten,250)

# Lower plot
p = [0.20, 0.20, 0.05, 0.52]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
ax0 = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
ax0.plot(traj.time[r_min:r_max]/1.e3, toten, zorder=1, **styles[0])
ax0.plot(traj.time[r_min:r_max]/1.e3, toten_ra, zorder=1, **styles[2])
#ax0.plot(traj.time[r_min:r_max]/1.e3, poten_ra, zorder=1, **styles[1])
ax0.set_xlabel('Time [ps]')
ax0.set_ylabel('$E$ [eV]')
#ax0.set_yticks(np.arange(-340,-300,10))
ymean = np.mean(toten[-500:-1])
ax0.set_ylim(ymean-10, ymean+10)

for line in ax0.yaxis.get_ticklines(): # line is a Line2D instance
    line.set_markersize(3)
for line in ax0.xaxis.get_ticklines(): # line is a Line2D instance
    line.set_markersize(3)


t_min = 500   # skip the first 500 equilibration steps
t_max = r_max

# Upper plot:

kB = 8.617343e-5 # Boltzmann constant in eV/K
deg_fr = 186 # degrees of freedom
temp = 2*ekin/deg_fr/kB

print "Temperature:"
temp_mean = np.mean(temp[t_min:t_max])
print "   Mean value:",temp_mean
temp_med = np.median(temp[t_min:t_max])
print "   Median value:",temp_med
temp_std = np.std(temp[t_min:t_max], dtype=np.float64)
print "   Standard deviation:",temp_std
temp_relvar = (temp_std**2)/(temp_med**2)
print "   Relative variance:",temp_relvar

temp_ra = plotutils.symmetric_running_mean(temp,250)

temp_mean_rounded = round(temp_mean/100.)*100
ax_temp_min = temp_mean_rounded - 550
ax_temp_max = temp_mean_rounded + 550
#ax_temp_min = 950 #temp_med - 5*temp_std
#ax_temp_max = 1950 #temp_med + 5*temp_std


p = [0.20, 0.59, 0.35, 0.05]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
ax1 = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
ax1.set_ylabel('$T$ [K]')
ax1.set_ylim(ax_temp_min, ax_temp_max)
ax1.plot(time[r_min:r_max]/1.e3, temp[r_min:r_max], zorder=1, **styles[0])
ax1.plot(time[r_min:r_max]/1.e3, temp_ra[r_min:r_max] , zorder=2, **styles[2])
ax1.fill_between(time[t_min:t_max]/1.e3, temp_med-temp_std, temp_med+temp_std, facecolor='yellow', alpha=0.5, zorder=0, linewidth=0)
for line in ax1.yaxis.get_ticklines(): # line is a Line2D instance
    line.set_markersize(3)
for line in ax1.xaxis.get_ticklines(): # line is a Line2D instance
    line.set_markersize(3)

p = [0.66, 0.59, 0.05, 0.05]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
ax2 = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
ax2.hist(temp[t_min:t_max], histtype='stepfilled',bins = 45, range = [ax_temp_min, ax_temp_max], facecolor='gray', edgecolor='gray', alpha=0.5, normed = True, linewidth=0.5, orientation='horizontal')
ax2.set_ylim(ax_temp_min, ax_temp_max)
ax2.set_xlim(0,0.008)
ax2.set_yticklabels([])
ax2.set_xticks([])
ax2.text(0.95,0.95,'$\sigma$ = %.f K' % (temp_std), horizontalalignment='right',
     verticalalignment='top',
     transform = ax2.transAxes
        )
for line in ax2.yaxis.get_ticklines(): # line is a Line2D instance
    line.set_markersize(3)


y=np.arange(ax_temp_min, ax_temp_max,1)
from matplotlib.mlab import normpdf
x=normpdf(y,temp_med,temp_std)
ax2.plot(x,y, color='black', linewidth = 1.)



#############################################################################
# (3) Save PDF 

fig_filename = os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.pdf'
sys.stdout.write("Writing %s... " % os.path.basename(fig_filename))
sys.stdout.flush()
plt.savefig(fig_filename)
sys.stdout.write("done!\n")

print

