#!/bin/bash
# lstart.sh (Mac, tmux 'yspec'): when yspec-l5, then yspec-p2 start running, start the production chunk in them
# (j has no autorun: a slot only runs what `j run` gives it). Restart this after a Mac reboot; it is idempotent per slot:
# a slot already given its run (marker ~/dns-yspec/mpi/yao/.started-SLOT) is skipped.
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; M=~/dns-yspec/mpi/yao
cd ~/science/xpred-collapse || exit 1
for s in "yspec-l5 $Y/l5.sh" "yspec-p2 $Y/prod4096.sh"; do
  set -- $s
  [ -e $M/.started-$1 ] && { echo "$(date +%T) $1 already started"; continue; }
  until ~/j/j up -w $1; do echo "$(date +%T) j up -w $1 failed, retrying"; sleep 60; done
  ~/j/j -n $1 run $2 && touch $M/.started-$1 && echo "$(date +%T) $1: started $2"
done
