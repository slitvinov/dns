#!/bin/bash
# psub.sh: Mac-side one-shot launcher of the approved Re 40000 4096 run: every 10 min, once the inputs exist (R1 r9kpad/done,
# R2 r9k768noS/done, G1 r20k2048p/done, a 4096 rehearsal report), runs src/gate.py on rc; on gate.ok submits three 5-node
# chunk slots yspec-p1..p3 (afterany chain, 6 h each, ban list excluded) running prod4096.sh, then exits.
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; PY=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/venv/bin/python
L=/private/tmp/claude-501/-Users-lisergey-web/b4e6832d-44e5-469b-bfc6-ba66a1c8ac03/scratchpad
cd ~/science/xpred-collapse || exit 1
until ssh rc "test -e $Y/gate.ok" 2>/dev/null; do
  ssh rc "cd $Y && test -e run/r9kpad/done && test -e run/r9k768noS/done && test -e run/r20k2048p/done && grep -l 'grid 4096' reh.*.txt > /dev/null 2>&1 && $PY src/gate.py > gate.log 2>&1" 2>/dev/null
  ssh rc "test -e $Y/gate.ok" 2>/dev/null && break
  sleep 600
done
echo "$(date +%T) gate.ok: $(ssh rc "cat $Y/gate.ok" 2>/dev/null | tr '\n' ' ')"
ex=$(grep -v '^#' ~/scienceclaw/gpu-banlist.txt | grep -o '^holygpu[0-9a-z]*' | sort -u | tr '\n' ',' | sed 's/,$//')
dep=""
for p in p1 p2 p3; do
  ~/j/j up yspec-$p -p gpu_h200,seas_gpu,gpu_requeue --constraint=h200 --exclude=$ex $dep -N 5 --ntasks-per-node=4 --cpus-per-task=8 --gres=gpu:4 --mem=400G -t 6:00:00 --signal=B:USR1@600
  id=$(ssh rc "squeue -u \$USER -h -o '%i %j' | grep 'yspec-$p\$' | cut -d' ' -f1" 2>/dev/null); dep="--dependency=afterany:$id"
  (~/j/j up -w yspec-$p && ~/j/j -n yspec-$p run $Y/prod4096.sh) > $L/w-$p.log 2>&1 &
  echo "$(date +'%F %H:%M') yspec-$p submitted (job $id, 5 nodes x4 H200, 6h, excl $ex): Re 40000 4096 production chunk prod4096.sh. Spectral Yao agent." >> OPS.md
done
wait
