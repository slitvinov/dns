#!/bin/bash
# prod4096.sh: one chunk of the Yao Re 40000 production run (Sergey's approval 2026-10-31): 4096^3 on 20 H200 (5 nodes),
# t* 1.3 -> 2.3 (t 2.7911 -> 4.9381), from the 1536 base dump at t* 1.3 (run/l40kbase, zero-padded), then from the newest
# k. dump of run/y40k4096. Binary b6, -S, C 2, -F 0.00015625, NCCL 8/32 channels, spectra every 0.02 in t, diagnostics
# also every 0.002 in t* over 1.75-1.95, snapshots (|w|, w_z, lambda2/1000, quarter crop x < 1.5, ~49 GB each; moved to the
# Mac and claw by mover.sh) at t* 1.5 1.6 1.7 1.75 1.8 1.85 1.9 2.0 2.15 2.3, k. dumps (413 GB) every 0.15 in t: a dump is
# deleted only after a newer one is complete (at most two on disk). First-steps guard: expected 2.3 s/step; a run slower
# than 3x is killed, its nodes go to badnodes.txt and are excluded from the pending yspec-p* chunks. Requires yspec/gate.ok.
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; D=$Y/run/y40k4096
[ -e $Y/gate.ok ] && [ ! -e $Y/hold ] || { echo "prod4096: no gate.ok"; scancel $SLURM_JOB_ID; exit 1; }
[ -e $D/finished ] && { echo "prod4096: finished"; scancel $SLURM_JOB_ID; exit 0; }
mkdir -p $D; [ -e $D/start ] || ln -s $(ls $Y/run/l40kbase/k.[0-9]*[0-9] | sort | tail -1) $D/start
KT=$(awk 'BEGIN{s = ""; for (t = 2.95; t < 4.93; t += 0.15) s = s (s == "" ? "" : ",") t; print s}')
W=3.2205,3.4352,3.6499,3.7572,3.8646,3.9719,4.0793,4.2939,4.6160,4.9381
( while sleep 300; do
    ls $D/k.[0-9]*[0-9] 2>/dev/null | while read f; do [ -L $f ] || echo $f; done | sort | head -n -1 | while read f; do rm -f $f $f.t; done
  done ) &
w=$!
before=$(cat $Y/badnodes.txt 2>/dev/null | wc -l)
EXPSPS=2.3 BIN=$Y/src/boxsymmp_b6 $Y/seg2.sh $D 4096 20 40000 2.0 4.9381 $KT -S -F 0.00015625 -e 0.02 -D 0.004294,3.75720,4.18660 -W $W
s=$?; kill $w
if [ $(cat $Y/badnodes.txt 2>/dev/null | wc -l) -gt $before ]; then
  ex=$( (echo holygpu8a12204; echo holygpu8a12402; tail -1 $Y/badnodes.txt | grep -o "holygpu[0-9a-z]*") | sort -u | tr "\n" "," | sed "s/,$//")
  for j in $(squeue -u $USER -h -t PD -o "%i %j" | grep "yspec-p" | cut -d" " -f1); do scontrol update jobid=$j ExcNodeList=$ex; done
  echo "$(date +%F\ %T) prod4096: slow nodes, pending chunks now exclude $ex" >> $D/seg.log
fi
[ $s = 0 ] && touch $D/finished
scancel $SLURM_JOB_ID
