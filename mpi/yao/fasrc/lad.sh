#!/bin/bash
# lad.sh N NP: G2 rung: Re 40000 on an N^3 grid (NP ranks) from the 1536 k. dump at t* 1.3 (run/l40kbase, zero-padded)
# to t 4.1866 (t* 1.95), C 2, -F 0.00125*512/N, -S, spectra every 0.05; no dumps. Waits for the base dump.
. /etc/profile > /dev/null 2>&1
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; . $Y/src/mpenv.sh
N=$1; NP=$2; B=$Y/src/boxsymmp_b4
until [ -e $Y/run/l40kbase/done ]; do sleep 60; done
k=$(ls $Y/run/l40kbase/k.[0-9]*[0-9] | sort | tail -1); s=$(basename $k)
D=$Y/run/l40k$N; mkdir -p $D && cd $D || exit 1
ln -sf $k $s; ln -sf $k.t $s.t
F=$(awk -v n=$N "BEGIN{printf \"%.6g\", 0.00125 * 512 / n}")
echo "$(date +%F\ %T) lad.sh n $N np $NP from $s F $F host $(hostname) $(scontrol show hostnames $SLURM_JOB_NODELIST 2>/dev/null | tr "\n" " ")" | tee -a lad.log
mpirun -np $NP --oversubscribe --bind-to none $B -x $N -y $N -z $N -X 2 -Y 2 -Z 2 -n 1.920096311e-4 -s 0.01 -c 2 -F $F -f -S \
  -e 0.02 -D 0.004294,3.75720,4.18660 -r $s -G 1536,1536,1536 -t 4.1866 > out 2> err
echo "$(date +%F\ %T) lad.sh exit $? $(grep -h s/step err)" | tee -a lad.log
