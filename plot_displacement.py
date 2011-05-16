#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import os, sys
import numpy as np
import getopt
from oppvasp.vasp.parsers import read_trajectory
from oppvasp.plotutils import DisplacementPlot

################################# OPTIONS ########################################

def usage():
    print "Usage: plot_displacement.py [first_step] [last_step] [atom] [filename]\n" \
            + "Generates a plot of r^2(t)"


try:    #parse the arguments
    opts, args = getopt.getopt(sys.argv[1:], "h", ["help", "first_step=","last_step=","atom=","filename="])
except getopt.GetoptError:
    usage()
    sys.exit(2)

first_step = 0
last_step = -1
atom_no = 0
filename = os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.pdf'

if len(args) > 0:
    first_step = int(args[0])
    if len(args) > 1:
        last_step = int(args[1])
        if len(args) > 2:
            atom_no = int(args[2])
            if len(args) > 3:
                filename = args[3]

for opt in opts:
    if opt[0] == '--filename':
        filename = opt[1]
    if opt[0] == '--first_step':
        first_step = int(opt[1])
    if opt[0] == '--last_step':
        last_step = int(opt[1])
    if opt[0] == '--atom':
        atom_no = int(opt[1])

print "Range: %i:%i. Atom: %i" % (first_step, last_step, atom_no)


############################## MAIN SCRIPT ######################################

# Read trajectory:
traj = read_trajectory(unwrap_pbcs = True, poscar_file = 'POSCAR')
traj.set_selection(first_step, last_step)

# Make plot
dp = DisplacementPlot(traj)
dp.add_plot( what = 'r2', smoothen = True, style = { 'linestyle' : '--', 'dashes' : (3,1), 'color' : 'black' } ) # avg r^2 for all atoms
#dp.add_plot( what = 'x', atom_no = 0, smoothen = True, linear_fit = True)
dp.add_plot( what = 'r2', atom_no = 0, smoothen = True, linear_fit = False)
dp.add_plot( what = 'r2', atom_no = 0, smoothen = True, linear_fit = True)
#dp.add_plot( what = 'x', atom_no = 0, smoothen = False, style = { 'zorder' : -1, 'alpha': 0.4, 'color': 'gray' } )

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
dp.ax1.set_ylim(ymin,ymax)

dp.save_plot(filename)

