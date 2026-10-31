#!/bin/bash
# g5.sh V1 [V2]: G1-cause tests at Re 20000, 1536^3 on 4 H200, from the 1536 k. dump at t* 1.3 (run/r20kbase) to t* 2.0
# (t 4.2939), diagnostics every 0.002 in t* over 1.75-1.95. Variants (reference = r20k1536: C 2, -S, -F 0.000416667):
#  nofilt: 2/3-rule mask (no -f), C 2, -S;  c1F: C 1, -S, -F;  c2pf: C 2, -S, filter once per step (no -F);
#  c1val: C 1, -S, filter once per step (the flags of the validated Re 5000/9000 runs)
. /etc/profile > /dev/null 2>&1
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; . $Y/src/mpenv.sh; B=$Y/src/boxsymmp_b5
for v in "$@"; do
  case $v in nofilt) X="-c 2 -S";; c1F) X="-c 1 -S -f -F 0.000416667";; c2pf) X="-c 2 -S -f";; c1val) X="-c 1 -S -f";; esac
  D=$Y/run/r20kt_$v; mkdir -p $D && cd $D || exit 1
  ln -sf $Y/run/r20kbase/k.00003532 k.00003532; ln -sf $Y/run/r20kbase/k.00003532.t k.00003532.t
  mpirun -np 4 --oversubscribe --bind-to none $B -x 1536 -y 1536 -z 1536 -X 2 -Y 2 -Z 2 -n 3.840192622e-4 -s 0.01 $X \
    -e 0.05 -D 0.004294,3.75720,4.18660 -r k.00003532 -t 4.2939 > out 2> err
  echo "g5 $v exit $? $(grep -h s/step err)" > done
done
scancel $SLURM_JOB_ID
