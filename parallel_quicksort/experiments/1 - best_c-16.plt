set terminal windows font "Tahoma,23"
set key nobox right top font ",22"
set xlabel 'c (cutoff = n/2^c)' offset 0,1
set ylabel 'speedup' offset 4.5,0
set xtics 1,1,23 font ",19" offset 0,0.5
set ytics font ",19" offset 0.7,0
set size ratio 0.5
set grid

set linetype 1 linewidth 3 linecolor rgb "#8A2BE2"       pointtype 12
set linetype 2 linewidth 3 linecolor rgb "red"           pointtype 9
set linetype 3 linewidth 3 linecolor rgb "blue"          pointtype 11
set linetype 4 linewidth 3 linecolor rgb "coral"         pointtype 6
set linetype 5 linewidth 3 linecolor rgb "magenta"       pointtype 2
set linetype 6 linewidth 3 linecolor rgb "forest-green"  pointtype 5
set linetype 7 linewidth 3 linecolor rgb "#6495ED"       pointtype 7   # cornflowerblue
set linetype 8 linewidth 3 linecolor rgb "olive"         pointtype 13
set linetype 9 linewidth 3 linecolor rgb "navy"          pointtype 4

set label 'm�ximo' at 3.3,1.9 font ',17'
set arrow from 4,2 to 4,2.4 linewidth 3
set label 'melhor c' at 3.3,0.6 font ',17'
set arrow from 4,0.5 to 4,0.1 linewidth 3


plot 'c16.dat' using 1:4 t 'n=2^{16}' with linespoints ls 2
     
     