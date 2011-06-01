
reset
set terminal X11 enhanced
set output
set xlabel "E - E_{Fermi} (eV)"
set ylabel "DOS (states/eV)"

plot 'DOS0' using 1:2 with lines title "Total DOS"
pause -1