#!/bin/bash
#SBATCH -J yspec-wstart -p shared -c 1 --mem=1G -t 2-00:00:00 -o /n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec/wstart.log
# wstart.sh: cluster-side starter for the y40k4096 chunks (no Mac needed): when the j slot yspec-l5 (then yspec-p2) is
# RUNNING, start its run through the slot's j daemon, as `j -n SLOT run` from the Mac would. prod4096.sh holds a per-job
# lock, so a second start (Mac lstart.sh too) exits at once. Ends when both are started or gone.
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec
for v in $(env | grep -o "^SLURM_[A-Z_]*"); do unset $v; done
go() { # go SLOT CMD: wait for RUNNING, then start CMD in the slot (once)
  local j st h
  while :; do
    read j st h <<< "$(squeue -u $USER -h -n j-$1 -o '%i %T %B' | head -1)"
    [ -z "$j" ] && { echo "$(date +%F\ %T) $1: no such job, skipped"; return; }
    [ "$st" = RUNNING ] && break; sleep 60
  done
  ls -d $Y/run/lock.$j.* > /dev/null 2>&1 && { echo "$(date +%F\ %T) $1 ($j): already started"; return; }
  sleep 60   # the j daemon's socket
  for i in 1 2 3 4 5; do
    srun --jobid=$j --overlap -n1 -N1 -w $h -Q ~/jcmd/j -sock /tmp/j-$(id -u)-$1.sock run $2 && { echo "$(date +%F\ %T) $1 ($j on $h): started $2"; return; }
    sleep 60
  done
  echo "$(date +%F\ %T) $1 ($j): FAILED to start $2"
}
go yspec-l5 $Y/l5.sh
go yspec-p2 $Y/prod4096.sh
