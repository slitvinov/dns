#!/bin/bash
# reh.sh NP: dress rehearsal of the Re 40000 fine phase on NP ranks (16 = 4 nodes), production config: binary b5 (= b4 +
# BOXSYM_RB, default off), -S, C 2, -F 0.000208333, NCCL 8/32 channels, from the t* 1.3 base dump padded to 3072.
#  A: ~40 steps with a k. dump after ~30 (-T), timed;  B: a new mpirun restarts from that dump to the same end time:
#  its out rows must equal A's;  C: exchanges alone, compute alone;  D: the slow-node guard with expected 0.05 s/step.
# Writes yspec/reh.<jobid>.txt, deletes its dump, cancels the slot.
. /etc/profile > /dev/null 2>&1
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; . $Y/src/mpenv.sh; B=$Y/src/boxsymmp_b5; NP=${1:-16}
O=$Y/reh.$SLURM_JOB_ID.txt; D=$Y/run/reh3072; rm -rf $D; mkdir -p $D && cd $D || exit 1
NC="-x NCCL_NCHANNELS_PER_NET_PEER=8 -x NCCL_MIN_NCHANNELS=32"
mp() { mpirun -np $NP --bind-to none $NC "$@"; }
k=$(ls $Y/run/l40kbase/k.[0-9]*[0-9] | sort | tail -1); s=$(basename $k); ln -s $k $s; ln -s $k.t $s.t
G="-x 3072 -y 3072 -z 3072 -X 2 -Y 2 -Z 2 -n 1.920096311e-4 -s 0.01 -c 2 -F 0.000208333 -f -S"
echo "$(date) nodes $(scontrol show hostnames $SLURM_JOB_NODELIST | tr '\n' ' ') np $NP binary b5" > $O
# A: from t 2.7911; dt at 3072 ~ 3.5e-4 -> dump at +0.0105 (~30 steps), end at +0.014 (~40 steps)
t1=2.8016; t2=2.8051
mp $B $G -k -e 1 -r $s -G 1536,1536,1536 -T $t1 -t $t2 > outA 2> errA &
p=$!; $Y/src/spscheck.sh outA 2.0 $p & wait $p; echo "A exit $?: $(grep -h 's/step\|wrote\|read' errA | tr '\n' ' ')" >> $O
grep -h "rank 0/\|rank 15/" errA >> $O
kd=$(ls k.[0-9]*[0-9] | grep -v $s | sort | tail -1); echo "dump $kd size $(stat -c %s $kd) B (estimate 3 x 1537^3 x 16 = $((3 * 1537 * 1537 * 1537 * 16)))" >> $O
# B: fresh mpirun from the dump to t2
mp $B $G -k -e 1 -r $kd -t $t2 > outB 2> errB; echo "B exit $?: $(grep -h 's/step\|read' errB | tr '\n' ' ')" >> $O
sa=$(basename $kd | sed 's/k\.0*//')
awk 'NR == FNR {if (FNR > 1) a[$1] = $0; next} FNR > 1 && $1 > s && ($1 in a) {n++; if (a[$1] != $0) d++} END {printf "restart: %d common out rows after step %s, %d differ\n", n, s, d + 0}' s=$sa outA outB >> $O
rm -f k.[0-9]*[0-9] k.[0-9]*[0-9].t 2>/dev/null; rm -f $s $s.t
# C: split, on the same nodes (bench.sh: -I yao, fixed dt)
BIN=$B $Y/bench.sh 3072 $NP 1 R_a2a BOXSYM_A2AONLY=5 -- -S; BIN=$B $Y/bench.sh 3072 $NP 10 R_comp BOXSYM_NOCOMM=1 -- -S
BIN=$B $Y/bench.sh 3072 $NP 20 R_full -- -S
grep -h "^R_.*np $NP " $Y/bench.txt | tail -3 >> $O
# D: guard
mp $B $G -r $Y/run/l40kbase/$s -G 1536,1536,1536 -t 3.5 > outD 2> errD &
p=$!; SPSLOG=$D/guard.txt $Y/src/spscheck.sh outD 0.05 $p; echo "guard exit $? (1 = killed): $(cat guard.txt 2>/dev/null)" >> $O
wait $p 2>/dev/null; echo "$(date) done" >> $O
scancel $SLURM_JOB_ID
