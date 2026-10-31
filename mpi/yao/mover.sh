#!/bin/bash
# mover.sh: moves finished Re 40000 snapshots (w.STEP_*.npy + w.STEP.txt) from FASRC run/y40k3072 to claw
# /data/share/render/yao/re40000/ through this Mac, checks the sizes, deletes them on netscratch. Polls every 10 min;
# stops when run/y40k3072/stop.mover exists.
D=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec/run/y40k3072; C=/data/share/render/yao/re40000
ssh claw "sudo mkdir -p $C" || exit 1
while :; do
  ssh rc "test -e $D/stop.mover" && break
  for s in $(ssh rc "cd $D 2>/dev/null && grep -ho 'snapshot w.[0-9]*' err.* | sed 's/snapshot //' | sort -u"); do
    for f in ${s}_wabs.npy ${s}_wz.npy ${s}_lam2.npy $s.txt; do
      n=$(ssh rc "stat -c %s $D/$f 2>/dev/null") || continue
      [ -z "$n" ] && continue
      ssh rc "cat $D/$f" | ssh claw "sudo tee $C/$f.part > /dev/null && sudo mv $C/$f.part $C/$f"
      m=$(ssh claw "sudo stat -c %s $C/$f 2>/dev/null")
      [ "$m" = "$n" ] && ssh rc "rm -f $D/$f" && echo "$(date +%T) moved $f ($n B)"
    done
  done
  V=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec/run/r9kvid; C9=/data/share/render/yao/re9000
  if ssh rc "test -e $V/done" && [ -z "$(ssh claw "ls $C9/README.txt 2>/dev/null")" ]; then
    ssh claw "sudo mkdir -p $C9"
    for f in $(ssh rc "cd $V/fields && ls *.npy"); do
      n=$(ssh rc "stat -c %s $V/fields/$f"); ssh rc "cat $V/fields/$f" | ssh claw "sudo tee $C9/$f.part > /dev/null && sudo mv $C9/$f.part $C9/$f"
      [ "$(ssh claw "sudo stat -c %s $C9/$f")" = "$n" ] && ssh rc "rm -f $V/fields/$f" && echo "$(date +%T) moved $f"
    done
    [ -z "$(ssh rc "ls $V/fields/*.npy 2>/dev/null")" ] && cat ~/dns-yspec/mpi/yao/README_re9000.txt | ssh claw "sudo tee $C9/README.txt > /dev/null" && echo "$(date +%T) re9000 fields done"
  fi
  ssh claw "sudo chown -R www-data:www-data /data/share/render" > /dev/null 2>&1
  sleep 600
done
