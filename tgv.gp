$grid << EOD
0064  #1f77b4  "tg, 64^3"
0128  #2ca02c  "tg, 128^3"
0256  #d62728  "tg, 256^3"
0512  #000000  "tg, 512^3"
EOD
$paper << EOD
0100  #999999  "not stated"
0200  #999999  "not stated"
0400  #2ca02c  "128^3"
0800  #2ca02c  "128^3 (inferred)"
1600  #d62728  "256^3"
3000  #d62728  "256^3"
EOD
ng = |$grid|
np = |$paper|
array N[ng]
array C[ng]
array L[ng]
array R[np]
array P[np]
array Q[np]
do for [i = 1:ng] {
    N[i] = word($grid[i], 1)
    C[i] = word($grid[i], 2)
    L[i] = word($grid[i], 3)
}
do for [j = 1:np] {
    R[j] = word($paper[j], 1)
    P[j] = word($paper[j], 2)
    Q[j] = word($paper[j], 3)
}
set term svg size 1200, 1200 font 'arial,20'
set output "img/tgv.svg"
set rmargin 12
set size sq
set xlabel "time"
set ylabel "rate of energy dissipation"
set ytics 0, 0.01, 0.02
$at << EOD
10.15  0.00528
10.15  0.00706
10.15  0.00935
10.15  0.01040
10.15  0.01127
10.15  0.01269
EOD
do for [j = 1:np] {
    set label j sprintf("Re = %d", R[j] + 0) at word($at[j], 1), word($at[j], 2) \
        textcolor rgb P[j] font 'arial,16'
    set label 10 + j sprintf("Re = %d: %s", R[j] + 0, Q[j]) at 0.3, 0.0195 - 0.0006 * j \
        textcolor rgb P[j] font 'arial,16'
}
set label 10 "grid of Brachet et al.:" at 0.3, 0.0195 font 'arial,16'
set key at graph 0.97, graph 0.97 right top
plot [0:10][0:0.02] \
     for [i = 1:ng] for [j = 1:np] sprintf("data/tg/%s/%s", N[i], R[j]) \
         u 2:(2 * column(4) / R[j]) w l lw 3 lc rgb C[i] t (j == 1 ? L[i] : ""), \
     for [j = 1:np] "img/ref.txt" u 1:($3 == R[j] && $4 == 1 ? $2 : 1/0) \
         w p pt 6 ps 1.2 lc rgb P[j] t (j == 1 ? "Brachet et al. (1983), figure 7" : ""), \
     for [j = 1:np] "img/ref.txt" u 1:($3 == R[j] && $4 == 0 ? $2 : 1/0) \
         w p pt 2 ps 1.0 lc rgb P[j] t (j == 1 ? "figure 7, curves cross (uncertain)" : "")
