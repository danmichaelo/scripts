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

vel = traj.get_velocities( coordinates = 'direct' )


#############################################################################
# (2) Plot

styles = [ { 'color' : c } for c in ['red','blue','black'] ]

prepare_canvas('10 cm') 
p = [0.15, 0.15, 0.05, 0.05]  # margins: left, bottom, right, top. Height: 1-.15-.46 = .39
fig = plt.figure()
ax1 = fig.add_axes([ p[0], p[1], 1-p[2]-p[0], 1-p[3]-p[1] ])
ax1.set_xlabel(r'Time [ps]')
#ax1.set_ylabel(ur'Si-P bond length [Å]')

#dp.ax1.legend(['%d' % (at) for at in include_atoms], loc='upper right', ncol=1, frameon = False)
#dp.ax1.set_ylabel(u'$r^2$ [Å]')

p_vel = vel[:,0]
tlen = np.min((10000,p_vel.shape[0]))

natoms = vel.shape[1]
print "Found %d atoms" % (natoms)

vacf0 = np.sum(np.sum(vel[1] * vel[1], axis=1)) * 1./natoms
auto_correl = np.zeros((tlen))
fac = 1./vacf0  # normalization factor
fac = fac * 1./natoms
for i in range(1,tlen):
    #auto_correl[i] = np.dot(p_vel[1], p_vel[i]) * fac
    vacf = np.sum(np.sum(vel[1] * vel[i], axis=1)) * fac
    auto_correl[i] = vacf
    #print i,auto_correl[i]


#acs = symmetric_running_mean(auto_correl,3000)
acs = auto_correl #disable smoothing

#ax1.plot(traj.time[:]/1.e3, vel[:,0,0]) # x cord of at 0
ax1.plot(traj.time[:tlen]/1.e3, acs) # velocity auto-correlation function


#plt.legend(lines, labels, ncol = 1, loc = 'upper left', frameon = False)
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




