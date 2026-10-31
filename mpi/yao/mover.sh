#!/bin/bash
# mover.sh: moves finished Yao fields off FASRC netscratch. Master copy on this Mac, then a copy on claw for Rendering
# (claw is temporary), then delete on netscratch, every step checked by size. Polls every 10 min; stops when
# yspec/stop.mover exists.
#  Re 40000 snapshots: rc run/y40k4096 (and y40k3072) w.STEP_{wabs,wz,lam2}.npy + w.STEP.txt -> ~/dns/bench/ref/yao/re40000/
#     -> claw /data/share/render/yao/re40000/   (only snapshots whose 'snapshot w.STEP' line is in err.*: complete)
#  Re 9000 fields: rc run/r9kvid/fields/*.npy (after run/r9kvid/done) -> ~/dns/bench/ref/yao/re9000/ -> claw .../re9000/
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec/run; M=~/dns/bench/ref/yao; C=/data/share/render/yao
mkdir -p $M/re40000/y40k4096 $M/re40000/y40k3072 $M/re9000 $M/spec/run
mv1() { # mv1 rcpath macdir clawdir
  local f=$(basename $1) n m k
  n=$(ssh rc "stat -c %s $1 2>/dev/null"); [ -z "$n" ] && return 1
  rsync -a --partial rc:$1 $2/ || return 1
  m=$(stat -f %z $2/$f); [ "$m" = "$n" ] || return 1
  ssh claw "sudo mkdir -p $3" && rsync -a --partial --rsync-path="sudo rsync" $2/$f claw:$3/ || return 1
  k=$(ssh claw "sudo stat -c %s $3/$f"); [ "$k" = "$n" ] || return 1
  ssh rc "rm -f $1" && echo "$(date +%T) moved $f ($n B): Mac + claw, deleted on netscratch"
}
while :; do
  ssh rc "test -e $Y/../stop.mover" && break
  for P in y40k4096 y40k3072; do
    for s in $(ssh rc "cd $Y/$P 2>/dev/null && grep -ho 'snapshot w.[0-9]*' err.* 2>/dev/null | sed 's/snapshot //' | sort -u"); do
      for f in ${s}_wabs.npy ${s}_wz.npy ${s}_lam2.npy $s.txt; do
        ssh rc "test -e $Y/$P/$f" && mv1 $Y/$P/$f $M/re40000/$P $C/re40000/$P
      done
    done
  done
  if ssh rc "test -e $Y/r9kvid/done"; then
    for f in $(ssh rc "cd $Y/r9kvid/fields 2>/dev/null && ls *.npy 2>/dev/null"); do mv1 $Y/r9kvid/fields/$f $M/re9000 $C/re9000; done
    if [ ! -e $M/re9000/README.txt ] && [ -z "$(ssh rc "ls $Y/r9kvid/fields/*.npy 2>/dev/null")" ] && ls $M/re9000/*.npy > /dev/null 2>&1; then
      cp ~/dns-yspec/mpi/yao/README_re9000.txt $M/re9000/README.txt
      rsync -a --rsync-path="sudo rsync" $M/re9000/README.txt claw:$C/re9000/ && echo "$(date +%T) re9000 fields done"
    fi
  fi
  # diagnostics of every run (out, err, e. spectra, logs, gate files) mirrored to the Mac (small; no fields, no dumps)
  rsync -a --prune-empty-dirs --include='*/' --include='out*' --include='err*' --include='e.*' --include='*.log' \
    --include='*.info' --include='*.t' --include='w.*.txt' --exclude='*' rc:$Y/ $M/spec/run/ > /dev/null 2>&1
  rsync -a rc:$Y/../'{gate.ok,gate.fail,gate.log,w4096.txt,w4096.go,w4096.no,gate4096.log,bench.txt}' $M/spec/ > /dev/null 2>&1
  ssh claw "sudo chown -R www-data:www-data /data/share/render" > /dev/null 2>&1
  sleep 600
done
