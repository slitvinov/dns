#!/bin/bash
# rvar.sh N V...: Re 9000 G1-cause tests on one GPU (~/yspec, stub MPI/NCCL build): base N^3 from t 0 to t* 1.3 with the
# validated flags (C 1, -S, Hou-Li once per step) and a k. dump; then variants from that dump to t* 2.0 (t 4.2939), with
# diagnostics every 0.002 in t* over 1.75-1.95. Variants: nofilt (2/3 rule, C 1, -S), c2pf (C 2, -S, -f), c1F (C 1, -S,
# -f -F 0.00125*512/N), c2F (C 2, -S, -f -F: the G1 / production flags), val (the validated flags, continued: reference).
cd ~/yspec || exit 1; N=$1; shift; B=$PWD/boxsymmp_b5; F=$(awk -v n=$N 'BEGIN{printf "%.6g", 0.00125 * 512 / n}')
G="-x $N -y $N -z $N -X 2 -Y 2 -Z 2 -n 8.533761e-4 -s 0.01"
mkdir -p r9b$N && cd r9b$N
[ -e done ] || { $B $G -c 1 -S -f -k -e 0.25 -I yao -T 2.7911 -t 2.7911 > out 2> err && touch done; }
k=$(ls k.[0-9]*[0-9] | sort | tail -1); cd ..
for v in "$@"; do
  case $v in nofilt) X="-c 1 -S";; c2pf) X="-c 2 -S -f";; c1F) X="-c 1 -S -f -F $F";; c2F) X="-c 2 -S -f -F $F";; val) X="-c 1 -S -f";; esac
  D=r9v${N}_$v; mkdir -p $D && cd $D; ln -sf ../r9b$N/$k $k; ln -sf ../r9b$N/$k.t $k.t
  $B $G $X -e 0.05 -D 0.004294,3.75720,4.18660 -r $k -t 4.2939 > out 2> err; echo "$v exit $? $(grep -h s/step err)" > done; cd ..
done
