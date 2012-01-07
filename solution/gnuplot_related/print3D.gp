set xlabel "x"
set ylabel "y"
set zlabel "z"
set xrange [0:2]
set yrange [0:2]
set zrange [0:2]
#pointtype 5: filled square
#pointtype 7: filled circle
#splot  'tmp3d.gpl'  using 1:2:3:4  with points pointtype 7 pointsize 2  palette
splot  'sol3d_refine4_Q2.gpl'  using 1:2:3:4  with points pointtype 7 pointsize 1  palette
pause -1
