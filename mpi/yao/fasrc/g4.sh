#!/bin/bash
# g4.sh: G1 2048 (Re 20000) quickly: 1536 from t 0 to t* 1.3 with a k. dump, zero-padded to 2048, to t* 2.0 (t 4.2939),
# C 2, -S, -F 0.00125*512/n, diagnostics also every 0.002 in t* over 1.75-1.95, spectra every 0.02 (run/r20k2048p)
. /etc/profile > /dev/null 2>&1
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; . $Y/src/mpenv.sh; B=$Y/src/boxsymmp_b4
D=$Y/run/r20k2048p; mkdir -p $D && cd $D || exit 1
mp() { mpirun -np 4 --oversubscribe --bind-to none "$@"; }
G="-X 2 -Y 2 -Z 2 -n 3.840192622e-4 -s 0.01 -f -c 2 -S"
[ -e k.done ] || { mp $B -x 1536 -y 1536 -z 1536 $G -F 0.000416667 -k -e 0.1 -I yao -T 2.7911 -t 2.7911 > out.0 2> err.0 && touch k.done; }
k=$(ls k.[0-9]*[0-9] | sort | tail -1)
mp $B -x 2048 -y 2048 -z 2048 $G -F 0.0003125 -e 0.02 -D 0.004294,3.75720,4.18660 -r $k -G 1536,1536,1536 -t 4.2939 > out.1 2> err.1
echo "g4.sh exit $? $(grep -h s/step err.1)" > done
rm -f k.[0-9]*
scancel $SLURM_JOB_ID
