# symbols: http://ayapin-film.sakura.ne.jp/Gnuplot/Docs/ps_guide.pdf

set terminal windows font "Tahoma,23"
set key nobox left top font ",22"
set xlabel 'itens' offset 0,1
set ylabel 'tempo (s)' offset 4.5,0
set format x '2^{%.0f}'
set xtics (16,17,18,19,20)  font ",19" offset 0,0.2
set ytics font ",19" offset 0.5,0
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


plot 'performance.dat' using 1:2 t 'SSqs'  with linespoints ls 1,\
     'performance.dat' using 1:9 t 'SOqs' with linespoints ls 4
     