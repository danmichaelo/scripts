#!/usr/bin/env python

import os,sys
import numpy as np

from optparse import OptionParser
from oppvasp.vasp.parsers import PoscarParser

parser = OptionParser()
parser.add_option('-d', '--doscar', dest = 'doscar', default = 'DOSCAR')
parser.add_option('-o', '--outcar', dest = 'outcar', default = 'OUTCAR')
(options, args) = parser.parse_args()

if not os.path.exists(options.doscar):
    print "Error: File %s not found" % (options.doscar)
    sys.exit(1)

if os.path.exists(options.outcar):
    fermi = [l for l in open(options.outcar) if 'E-fermi' in l ]
    fermi = float(fermi[-1].split()[2])
    print "Fermi level: %.5f eV " % (fermi)
else:
    print "Warning: File %s not found. The DOS will not be adjusted to the Fermi level" % (options.outcar)
    fermi = 0.0

doscar = open(options.doscar)
main_header = []
for i in range(5):
    main_header.append(doscar.readline())

natoms = int(main_header[0].split()[0])
#print "%d atoms" % (natoms)

n = 0
while True:
    try:
        # read section header line
        (emax,emin,npoints,unknown1,unknown2) = map(float,doscar.readline().split())
    except:
        # not a header line, we've reached EOF or some obscure line
        break
    npoints = int(npoints)
    firstline = doscar.readline().split()
    ncols = len(firstline)
    dos = np.zeros((npoints,ncols))
    dos[0] = firstline
    for i in range(1,npoints):
        dos[i] = doscar.readline().split()

    # Subtract Fermi level from the energy-column (the first column)
    dos[:,0] -= fermi

    np.savetxt('DOS_%d'%(n), dos)
    if n == 0:
        print "Wrote total DOS to DOS_%d (%d points from %.1f to %.1f)" % (n,npoints,emin,emax)
    n += 1
if n > 1:
    print "Wrote site-projected DOS to files DOS_[1-%d]" % (n-1)


