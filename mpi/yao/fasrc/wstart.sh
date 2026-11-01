#!/bin/bash
#SBATCH -J yspec-wstart -p sapphire,shared -c 1 --mem=1G -t 2-00:00:00 -o /n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec/wstart.log
# wstart.sh: cluster-side starter for the y40k4096 chunks (no Mac needed). For the j slots yspec-l5 and yspec-p2: each time a
# slot's job is RUNNING with a (job id, restart count) not yet started, start prod4096.sh in it through the slot's j daemon,
# as `j -n SLOT run` from the Mac would (so a requeued job is restarted too; prod4096 resumes from the newest verified dump).
# prod4096.sh holds lock.JOBID.RESTARTCOUNT, so a second start (Mac lstart.sh too) exits at once. Ends when the run is
# finished/stopped or both jobs are gone.
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; D=$Y/run/y40k4096; L=$D/seg.log
for v in $(env | grep -o "^SLURM_[A-Z_]*"); do unset $v; done
log() { echo "$(date +%F\ %T) wstart: $*"; mkdir -p $D; echo "$(date +%F\ %T) wstart: $*" >> $L; }
log "watching yspec-l5, yspec-p2"
while :; do
  [ -e $D/finished ] && { log "run finished, exit"; exit 0; }
  [ -e $D/stopped ] && { log "run stopped, exit"; exit 0; }
  n=0
  for s in yspec-l5 yspec-p2; do
    read j st h <<< "$(squeue -u $USER -h -n j-$s -o '%i %T %B' | head -1)"
    [ -z "$j" ] && continue
    n=$((n + 1)); [ "$st" = RUNNING ] || continue
    r=$(scontrol show job $j | grep -o "Restarts=[0-9]*" | cut -d= -f2)
    [ -d $Y/run/lock.$j.${r:-0} ] && continue
    [ -e /tmp/wstart.$j.${r:-0} ] && continue        # already tried this one
    sleep 60                                          # the j daemon's socket
    ok=0
    for i in 1 2 3 4 5; do
      srun --jobid=$j --overlap -n1 -N1 -w $h -Q ~/jcmd/j -sock /tmp/j-$(id -u)-$s.sock run $Y/prod4096.sh && { ok=1; break; }
      sleep 60
    done
    touch /tmp/wstart.$j.${r:-0}
    [ $ok = 1 ] && log "$s job $j restart ${r:-0} on $h: started prod4096.sh" || log "$s job $j restart ${r:-0}: FAILED to start prod4096.sh"
  done
  [ $n = 0 ] && { log "no yspec-l5/p2 jobs left, exit"; exit 0; }
  sleep 60
done
