#!/bin/bash
# prod4096.sh: one chunk of the Yao Re 40000 production run (Sergey's approval 2026-10-31): 4096^3 on 20 H200 (5 nodes),
# t* 1.3 -> 2.3 (t 2.7911 -> 4.9381), from the 1536 base dump at t* 1.3 (run/l40kbase, zero-padded), then from the newest
# k. dump of run/y40k4096. Binary b7, -S, C 2, -F 0.00015625, NCCL 8/32 channels, spectra every 0.02 in t, diagnostics
# also every 0.002 in t* over 1.75-1.95, snapshots (|w|, w_z, lambda2/1000, quarter crop x < 1.75, ~57 GB, plus the whole quarter block-averaged 4x4 in x, y
# (d.STEP_*, 2 x 2.2 GB); moved to the
# Mac and claw by mover.sh) at t* 1.5 1.6 1.7 1.75 1.8 1.85 1.9 2.0 2.15 2.3, k. dumps (413 GB) every 0.15 in t: a dump is
# deleted only after a newer one is complete (at most two on disk). First-steps guard: expected 2.3 s/step; a run slower
# than 3x is killed, its nodes go to badnodes.txt and are excluded from the pending yspec-p* chunks. Requires yspec/gate.ok.
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; D=$Y/run/y40k4096
[ -e $Y/gate.ok ] && [ ! -e $Y/hold ] || { echo "prod4096: no gate.ok"; scancel $SLURM_JOB_ID; exit 1; }
[ -e $D/finished ] && { echo "prod4096: finished"; scancel $SLURM_JOB_ID; exit 0; }
[ -e $D/stopped ] && { echo "prod4096: stopped (see seg.log)"; scancel $SLURM_JOB_ID; exit 1; }
mkdir -p $D; [ -e $D/start ] || ln -s $(ls $Y/run/l40kbase/k.[0-9]*[0-9] | sort | tail -1) $D/start
# k. dumps: an early one at t 2.85 (verifies the write in the first chunk), then every 0.15 in t
KT=$(awk 'BEGIN{s = "2.85"; for (t = 2.95; t < 4.93; t += 0.15) s = s "," t; print s}')
W=3.2205,3.4352,3.6499,3.7572,3.8646,3.9719,4.0793,4.2939,4.6160,4.9381
# an older dump is deleted only when the newest one is complete: renamed from .tmp (after fsync), full size, with its .t
( while sleep 120; do
    # 900 GB guard: while a dump is being written, if yspec is above 900 GB drop the older complete dump now
    if ls $D/k.*.tmp > /dev/null 2>&1 && [ "$(du -s --block-size=1G $Y | cut -f1)" -gt 900 ]; then
      o=$(ls $D/k.[0-9]*[0-9] 2>/dev/null | while read f; do [ -L $f ] || echo $f; done | sort | head -1)
      [ -n "$o" ] && rm -f $o $o.t && echo "$(date +%F\ %T) prod4096: yspec > 900 GB during a dump write: deleted $o early" >> $D/seg.log
    fi
    n=$(ls $D/k.[0-9]*[0-9] 2>/dev/null | while read f; do [ -L $f ] || echo $f; done | sort | tail -1)
    [ -n "$n" ] && [ -s $n.t ] && [ "$(stat -c %s $n)" = "$((3 * 2049 * 2049 * 2049 * 16))" ] || continue
    ls $D/k.[0-9]*[0-9] 2>/dev/null | while read f; do [ -L $f ] || echo $f; done | sort | head -n -1 | while read f; do rm -f $f $f.t; done
  done ) &
w=$!
before=$(cat $Y/badnodes.txt 2>/dev/null | wc -l)
# first-chunk check (Sergey): after ~20 min of stepping the s/step must be <= 1.5 x the 2.3 s model, else stop and report
( o=""; while [ -z "$o" ]; do sleep 30; o=$(ls -t $D/out.* 2>/dev/null | head -1); done
  st() { awk 'NR > 1 {s = $1} END {print s + 0}' $o; }; until [ "$(st)" -gt 0 ]; do sleep 20; done
  e=${o/out./err.}; wt() { grep -h "wrote '\|snapshot w" $e | sed "s/.* in \([0-9.]*\) s.*/\1/" | awk '{s += $1} END {print s + 0}'; }
  s0=$(st); t0=$(date +%s); w0=$(wt); sleep 1200; s1=$(st); t1=$(date +%s); w1=$(wt)
  sps=$(awk -v d=$((s1 - s0)) -v t=$((t1 - t0)) -v w0=$w0 -v w1=$w1 'BEGIN{printf "%.3f", d > 0 ? (t - (w1 - w0)) / d : 1e9}')
  echo "$(date +%F\ %T) prod4096 chunk check: $sps s/step over $((s1 - s0)) steps (limit 3.45); $(grep -h 'rank 0/' $D/err.* | tail -1 | grep -o 'gpu used [0-9.]* GB')" >> $D/seg.log
  awk -v s=$sps 'BEGIN{exit !(s > 3.45)}' && { echo "$(date +%F\ %T) prod4096: too slow, stopped" >> $D/seg.log; touch $D/stopped; pkill -u $USER -f "[b]oxsymmp_b7 -x 4096"; } ) &
cw=$!
EXPSPS=2.3 BIN=$Y/src/boxsymmp_b7 $Y/seg3.sh $D 4096 20 40000 2.0 4.9381 $KT -S -F 0.00015625 -e 0.02 -D 0.004294,3.75720,4.18660 -W $W -w 1.75 -V 4
s=$?; kill $w $cw 2> /dev/null
if [ $(cat $Y/badnodes.txt 2>/dev/null | wc -l) -gt $before ]; then
  ex=$( (echo holygpu8a12204; echo holygpu8a12402; tail -1 $Y/badnodes.txt | grep -o "holygpu[0-9a-z]*") | sort -u | tr "\n" "," | sed "s/,$//")
  for j in $(squeue -u $USER -h -t PD -o "%i %j" | grep "yspec-p" | cut -d" " -f1); do scontrol update jobid=$j ExcNodeList=$ex; done
  echo "$(date +%F\ %T) prod4096: slow nodes, pending chunks now exclude $ex" >> $D/seg.log
fi
[ $s = 0 ] && touch $D/finished
scancel $SLURM_JOB_ID
