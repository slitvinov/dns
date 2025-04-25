list="0400 0800 1600 3000"
set key off
plot [0:10][0:0.02] \
     for [r in list] \
     "0128/" . r u 2:(2*$4/r) w l lw 4 dt 7 lc 8, \
     for [r in list] \
     r u 2:(2*$4/r) w l lw 4 lc 8
