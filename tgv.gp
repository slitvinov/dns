$grid << EOD
0128  #2ca02c  7
0256  #d62728  5.5
0512  #000000  4
1024  #9467bd  2.5
EOD
$paper << EOD
0100  #999999
0200  #999999
0400  #2ca02c
0800  #2ca02c
1600  #d62728
3000  #d62728
EOD
ng = |$grid|
np = |$paper|
array N[ng]
array C[ng]
array W[ng]
array R[np]
array P[np]
do for [i = 1:ng] {
    N[i] = word($grid[i], 1)
    C[i] = word($grid[i], 2)
    W[i] = word($grid[i], 3)
}
do for [j = 1:np] {
    R[j] = word($paper[j], 1)
    P[j] = word($paper[j], 2)
}
set term svg size 1200, 1200 font 'arial,20'
set output "img/tgv.svg"
set key top left noenhanced
set size sq
set xlabel "time"
set ylabel "rate of energy dissipation"
set ytics 0, 0.01, 0.02
set xrange [0:10]
set yrange [0:0.02]
plot \
     for [i = 1:ng] for [j = 1:np] sprintf("data/tg/%s/%s", N[i], R[j]) \
         u 2:(2 * column(4) / R[j]) w l lw W[i] lc rgb C[i] \
         t (j == 1 ? sprintf("n = %d", N[i] + 0) : ""), \
     "data/tg/2048/3000" u 2:(2 * column(4) / 3000) w l lw 1.2 lc rgb "#ff7f0e" \
         t "n = 2048", \
     for [j = 1:np] "img/ref.txt" u 1:($3 == R[j] ? $2 : 1/0) \
         w p pt 6 lc rgb P[j] t (j == 1 ? "Brachet et al. (1983)" : "")
set output "img/tgv_zoom.svg"
set xrange [8:10]
set yrange [0.009:0.016]
set ytics 0.009, 0.001, 0.016
replot
set output "img/tgv_diff.svg"
set xrange [0:10]
set yrange [-0.0015:0.0015]
set ytics -0.0015, 0.0005, 0.0015
set ylabel "rate of energy dissipation minus n = 2048, Re = 3000"
diff(n) = sprintf("< awk 'NR == FNR {e[sprintf(\"%%.4f\", $2)] = $4; next} \
    (k = sprintf(\"%%.4f\", $2)) in e {print $2, 2 * ($4 - e[k]) / 3000}' \
    data/tg/2048/3000 data/tg/%s/3000", n)
plot for [i = 1:ng] diff(N[i]) w l lw W[i] lc rgb C[i] t sprintf("n = %d", N[i] + 0), \
     0 w l lw 1 lc rgb "#ff7f0e" t "n = 2048"
