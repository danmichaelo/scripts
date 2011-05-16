#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import os, sys
import numpy as np
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
     RotatingMarker, ReverseBar, SimpleProgress

from oppvasp.plotutils import prepare_canvas, symmetric_running_mean, symmetric_running_median
import matplotlib.pyplot as plt

from oppvasp.md import Trajectory
from oppvasp.vasp.parsers import IterativeVasprunParser

#############################################################################
# (1) Extract data


def read_trajectory( dir = './', xml_file ='vasprun.xml', npz_pbc_file = 'trajectory_pbc.npz', npz_file = 'trajectory_nopbc.npz', POSCAR_file = '' ):
    if os.path.isfile(dir + npz_file):
        traj = Trajectory(filename = dir +npz_file)
    else:
        p = IterativeVasprunParser(dir + xml_file)
        traj = p.get_all_trajectories()
        traj.save(dir + npz_pbc_file)
        if POSCAR_file != '':
            poscar = PoscarParser(dir + POSCAR_file)
            pos = poscar.get_positions( coordinates = 'direct' )
            print "Unwrapping using given initial pos"
            traj.unwrap_pbc( init_pos = pos)
        else:
            traj.unwrap_pbc()
        traj.save(dir + npz_file)
    return traj

traj = read_trajectory( dir = './' )

pos = traj.get_all_trajectories( coordinates = 'cartesian' )
pos = np.transpose(pos,[1,0,2]) # switch axes 0 and 1 for simple subtraction
r = pos - pos[0]
dist = np.sqrt(np.sum(r**2, axis=2))  # ** and sqrt are element-wise
dist = np.transpose(dist)  # switch axes back 
# dist.shape = (num_steps, num_atoms)

# Find nearest neighbours
arg = np.argsort(dist)
nearest = arg[0,1:3]

#sorted_dist = np.sort(dist) # shape: (num_steps, num_atoms)

#############################################################################
# (2) Plot

styles = [ { 'color' : c } for c in ['red','blue','black'] ]

prepare_canvas('10 cm') 
p = [0.15, 0.15, 0.05, 0.05]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
fig = plt.figure()
ax1 = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
ax1.set_xlabel(r'Time [ps]')
ax1.set_ylabel(ur'Si-P bond length [Å]')
ax1.set_ylim(2.0,3.0)
plotdata = []

#dp.ax1.legend(['%d' % (at) for at in include_atoms], loc='upper right', ncol=1, frameon = False)
#dp.ax1.set_ylabel(u'$r^2$ [Å]')

labels = []
lines = []
for i in nearest:
    d = dist[:,i]
    ds = symmetric_running_mean( d , 500 )
    #ax1.plot(traj.time[:], d, color = 'gray', alpha = 0.5)
    #line = ax1.plot(traj.time[:]/1.e3, ds )
    line = ax1.plot(traj.time[:]/1.e3, d )
    lines.append(line[0])
    print line
    print d.shape
    avg = np.sum(d)/d.shape[0]
    #ax1.axhline(y=avg, linestyle='dashed')
    print i,avg
    labels.append('No.%d' % (i))

print u"Shortest Si-P bond length at end of run:",dist[-1,nearest[0]],"Å for Si atom",nearest[0]
si = pos[nearest[0],-1]
p = pos[0,-1]
print "Si:",si
print "P:",p
print (si+p)/2./10.94

plt.legend(lines, labels, ncol = 3, loc = 'lower left', frameon = False)
#D = dp.get_diffusion_coeff()[0]
#dp.ax1.text(0.95,0.05,'$D$ = %.1f cm$^2$/s' % (D), horizontalalignment='right', verticalalignment='bottom', transform=ax1.transAxes)
#print "Diffusion coefficient:",D,"cm^2/s"




#############################################################################
# (3) Save PDF 

fig_filename = os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.pdf'
sys.stdout.write("Writing %s... " % os.path.basename(fig_filename))
sys.stdout.flush()
plt.savefig(fig_filename)
sys.stdout.write("done!\n")


print




