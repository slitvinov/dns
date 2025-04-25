list="0100 0200 0400 0800 1600 3000"
set term svg size 1200, 1200 font 'arial,20'
set output "img/tgv.svg"
set key off
set size sq
set xlabel "time"
set ylabel "rate of energy dissipation"
set ytics 0, 0.01, 0.02
plot [0:10][0:0.02] \
     "img/ref.txt" w p lc 8 pt 6, \
     for [r in list] "0256/" . r u 2:(2*$4/r) w l lw 3 lc 8, \
