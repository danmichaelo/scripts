#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- 
# vim:fenc=utf-8:et:sw=4:ts=4:sts=4:tw=0

import os, sys
import numpy as np
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
     RotatingMarker, ReverseBar, SimpleProgress

from oppvasp.plotutils import DisplacementPlot

dp = DisplacementPlot()
dp.read_trajectory( dir = './' )

# Find the diffusivity of all atoms:
m = 65

D = np.zeros(m)
pbar = ProgressBar(widgets=['Henter diff koeff...',Percentage(),Bar()], maxval = m).start()    
for at in range(m):
    pbar.update(at)
    D[at] = dp.get_diffusion_coeff(at)[0]
pbar.finish()

# Find the 3 atoms with the highest diffusivity:
arg = np.argsort(D)
include_atoms = arg[:-4:-1] # return the last three values (not including -4)


#pbar = ProgressBar(widgets=['Smoothen trajectories...',Percentage(),Bar()], maxval = m).start()    
styles = [ { 'color' : c } for c in ['red','blue','black'] ]
i = 0
for at in include_atoms:
    #pbar.update(at)
    dp.add_plot( what = 'r2', atom_no = at, smoothen = True, style = styles[i] )
    i += 1
#pbar.finish()

#dp.add_plot( what = 'r2', atom_no = 1, smoothen = True, style = { 'color': 'black' } )

dp.ax1.legend(['%d' % (at) for at in include_atoms], loc='upper right', ncol=1, frameon = False)
dp.ax1.set_ylabel(u'$r^2$ [Å]')


#dp.ax1.legend(('All atoms', 'P' % (p1[0],p1[1])), ncol = 2, loc = 'upper left', frameon = False)
#D = dp.get_diffusion_coeff()[0]
#dp.ax1.text(0.95,0.05,'$D$ = %.1f cm$^2$/s' % (D), horizontalalignment='right', verticalalignment='bottom', transform=ax1.transAxes)
#print "Diffusion coefficient:",D,"cm^2/s"


filename = os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.pdf'
dp.save_plot(filename)


