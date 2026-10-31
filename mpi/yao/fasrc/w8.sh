#!/bin/bash
# w8.sh: the 4096^3 window (Sergey pre-approved), only on w4096.go: 4096 on 32 H200 from the kept 3072 k. dump at t 3.75
# (t* 1.747), zero-padded, to t* 1.95 (t 4.1866), C 2, -F 0.00015625, -S, -D as the rungs, spectra every 0.02, snapshot
# at t* 1.85; no k. dumps (a requeue restarts the window from t 3.75). Frees the rolling dump of the finished production
# and the base first. Cancels its slot at the end.
. /etc/profile > /dev/null 2>&1
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; . $Y/src/mpenv.sh; PY=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/venv/bin/python
P=$Y/run/y40k3072; D=$Y/run/w40k4096
tend() { awk 'NR > 1 {t = $2} END {print t * 0.4657717}' $(ls $P/out.* | sort | tail -1); }
[ -e $Y/w4096.go ] || [ -e $Y/w4096.no ] || { awk -v t=$(tend) 'BEGIN{exit !(t >= 1.95)}' && $PY $Y/src/gate4096.py > $Y/gate4096.log 2>&1; }
[ -e $Y/w4096.go ] || { echo "w8.sh: no go ($(ls $Y/w4096.* 2>/dev/null))"; scancel $SLURM_JOB_ID; exit 0; }
k=$(for f in $P/k.[0-9]*[0-9]; do awk -v f=$f '{if ($1 > 3.749 && $1 < 3.751) print f}' $f.t; done | head -1)
[ -n "$k" ] || { echo "w8.sh: no t 3.75 dump"; scancel $SLURM_JOB_ID; exit 1; }
n=$(ls $P/k.[0-9]*[0-9] | sort | tail -1); awk '{exit !($1 > 4.938)}' $n.t && [ $n != $k ] && rm -f $n $n.t && rm -f $Y/run/l40kbase/k.*
mkdir -p $D && cd $D || exit 1
s=$(basename $k); ln -sf $k $s; ln -sf $k.t $s.t
echo "$(date +%F\ %T) w8.sh: 4096 from $s on 32 ranks $(scontrol show hostnames $SLURM_JOB_NODELIST | tr '\n' ' ')" | tee -a w8.log
mpirun -np 32 --bind-to none $Y/src/boxsymmp_b4 -x 4096 -y 4096 -z 4096 -X 2 -Y 2 -Z 2 -n 1.920096311e-4 -s 0.01 -c 2 \
  -F 0.00015625 -f -S -e 0.02 -D 0.004294,3.75720,4.18660 -W 3.9719 -r $s -G 3072,3072,3072 -t 4.1866 > out 2> err
echo "$(date +%F\ %T) w8.sh: exit $? $(grep -h s/step err)" | tee -a w8.log
scancel $SLURM_JOB_ID
