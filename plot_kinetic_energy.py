#!/usr/bin/env python

from oppvasp.plotutils import symmetric_running_median, prepare_canvas
prepare_canvas('10 cm')
from oppvasp.vasp.parsers import read_trajectory
import numpy as np
import matplotlib.pyplot as plt

t = read_trajectory()
ekin = t.kinetic_energy/t.num_atoms*1.e3
avgEkin = np.median(ekin)
plt.plot(t.time/1000., ekin, color='blue', alpha=0.2, label = 'Raw data')
plt.plot(t.time/1000., symmetric_running_median(ekin, 250), color='blue', alpha=1, label = '0.25 ps symmetric median')
plt.xlabel("Time [ps]")
plt.ylabel("Kinetic energy per atom [meV]")
plt.axhline(y=avgEkin, color='black', linestyle='--', label = 'Whole run median (%.0f meV)' % avgEkin)
plt.legend(frameon=False,loc='lower right')
plt.savefig('plot_kinetic_energy.pdf')

