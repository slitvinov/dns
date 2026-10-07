$grid << EOD
0064  #1f77b4
0128  #2ca02c
0256  #d62728
0512  #000000
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
array R[np]
array P[np]
do for [i = 1:ng] {
    N[i] = word($grid[i], 1)
    C[i] = word($grid[i], 2)
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
plot [0:10][0:0.02] \
     for [i = 1:ng] for [j = 1:np] sprintf("data/tg/%s/%s", N[i], R[j]) \
         u 2:(2 * column(4) / R[j]) w l lw 3 lc rgb C[i] \
         t (j == 1 ? sprintf("n = %d", N[i] + 0) : ""), \
     for [j = 1:np] "img/ref.txt" u 1:($3 == R[j] ? $2 : 1/0) \
         w p pt 6 lc rgb P[j] t (j == 1 ? "Brachet et al. (1983)" : "")
