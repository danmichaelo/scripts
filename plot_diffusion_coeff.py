#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import os, sys
import numpy as np
from oppvasp.vasp.parsers import read_trajectory
from oppvasp.plotutils import DisplacementPlot

first_step = 0
last_step = -1
if len(sys.argv) > 1:
    first_step= int(sys.argv[1])
    print "Starting from step %i" % (first_step)
    if len(sys.argv) > 2:
        last_step = int(sys.argv[2])
        print "Ending at step %i" % (last_step)

# Read trajectory:
traj = read_trajectory(unwrap_pbcs = True, poscar_file = 'POSCAR')
traj.set_selection(first_step, last_step)

# Make plot
dp = DisplacementPlot(traj)
dp.add_plot( what = 'D', atom_no = 0, smoothen = True, linear_fit = False)

#dp.ax1.legend(('All atoms', 'P'), ncol = 1, loc = 'upper left', frameon = False)
D = dp.get_diffusion_coeff([0])[0]


td = np.zeros((D.shape[0],2))
td[:,0] = dp.trajs[0].time
td[:,1] = D
np.savetxt('out.csv',td, delimiter='\t')

Ds = D[-1]
Dm = np.mean(D[:])
dp.ax1.axhline(Dm, color='red')
dp.ax1.text(0.95,0.95,'$D_P$ = %.2e cm$^2$/s' % (Dm), horizontalalignment='right', verticalalignment='top', transform=dp.ax1.transAxes)
print "Diffusion coefficient:",Ds,"cm^2/s"
print "Mean diffusion coefficient:",Dm,"cm^2/s"

ymax = 0
ymin = 0
for p in dp.plotdata:
    ymax_tmp = np.max(p['y'])
    if ymax_tmp > ymax:
        ymax = ymax_tmp
    ymin_tmp = np.min(p['y'])
    if ymin_tmp < ymin:
        ymin = ymin_tmp
ymax = (ymax)*1.1
ymin = (ymin)*1.1

dp.ax1.set_ylim(1.e-7,1.e-4)
dp.ax1.set_yscale('log')

filename = os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.pdf'
dp.save_plot(filename)

