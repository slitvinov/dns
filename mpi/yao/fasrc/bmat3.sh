#!/bin/bash
# bmat3.sh NP N...: timing of boxsymmp_b3 (row offsets once per tile row): newS_b3, new_b3, comp_b3, a2a-only float
R=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; NP=$1; shift; K=${K:-20}
export BIN=$R/src/boxsymmp_b3
for N in "$@"; do
  $R/bench.sh $N $NP $K newS_b3 -- -S
  $R/bench.sh $N $NP $K newShc2_b3 BOXSYM_HC=2 -- -S
  $R/bench.sh $N $NP $K new_b3
  $R/bench.sh $N $NP $K comp_b3 BOXSYM_NOCOMM=1
done
