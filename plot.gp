list="0100 0200 0400 0800 1600 3000"
set key off
plot [0:10][0:0.02] \
     for [r in list] \
     "0064/" . r u 2:(2*$4/r) w l lw 6 lc 1, \
     for [r in list] \
     r u 2:(2*$4/r) w l lw 1 lc 2
