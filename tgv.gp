list="0100 0200 0400 0800 1600 3000"
set term svg size 1200, 1200 font 'arial,20'
set output "img/tgv.svg"
set key top left
set size sq
set xlabel "time"
set ylabel "rate of energy dissipation"
set ytics 0, 0.01, 0.02
set label "Re = 100" at 0.2, 0.0082
set label "Re = 200" at 0.2, 0.0045
set label "Re = 400" at 0.2, 0.0026
set label "Re = 800, 1600, 3000" at 0.2, 0.0016
plot [0:10][0:0.02] \
     "img/ref.txt" w p lc 8 pt 6 t "Brachet et al. (1983)", \
     for [r in list] "0256/" . r u 2:(2*$4/r) w l lw 3 lc 8 t (r eq "0100" ? "fourier, 256^3" : ""), \
     for [r in list] "tg0256/" . r u 2:(2*$4/r) w l lw 3 lc 7 dt 2 t (r eq "0100" ? "tg -M 128" : "")
