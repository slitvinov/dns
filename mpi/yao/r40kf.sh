#!/bin/bash
# r40kf.sh: Re 40000 t* 1.3 512^3 fields from the l40kbase low-mode slab: rc -> hal (k2fields) -> Mac -> claw; deletes temporaries
# only after every copy is size-checked
set -e
M=~/dns/bench/ref/yao/re40000; Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; H=/scratch/lisergey/yspec; C=/data/share/render/yao/re40000
P=re40000_tstar1.30
rsync -a rc:$Y/l40kbase_m512.slab $M/
rsync -a $M/l40kbase_m512.slab hal:$H/
ssh hal "cd $H && mkdir -p r40kf && OMP_NUM_THREADS=16 python3 k2fields.py l40kbase_m512.slab 1536 512 r40kf/$P > r40kf/log 2>&1; cat r40kf/log"
rsync -a "hal:$H/r40kf/${P}_*.npy" $M/
cat > $M/README_tstar1.30.txt <<'R'
re40000_tstar1.30_{wabs,wz,lam2}.npy: Yao & Hussain antiparallel tubes, Re = Gamma0/nu = 40000, t* = 1.30 (t = 2.7911),
from the boxsymmp 1536^3 run (run/l40kbase, k.00003549), spectrally truncated to 512^3 (modes |k| < 256), full 2 pi box.
Pre-peak and under-resolved at 1536 (dry-run data for the video matcher, not a converged field).
float16, shape (512, 512, 512), index order [i][j][k] = (x, y, z), x_i = -pi + i 2pi/512 (same for y, z):
x collision normal, y propagation, z tube axis. wabs = |omega|, wz = omega_z, lam2 = lambda2/1000 (lambda2 of S^2 + W^2,
negative in vortex cores). Made by k2fields.py (~/dns-yspec/mpi/yao). max|w| 32.75 at 512^3; min lambda2 -172.3.
R
ssh claw "sudo mkdir -p $C"
rsync -a --rsync-path="sudo rsync" $M/${P}_wabs.npy $M/${P}_wz.npy $M/${P}_lam2.npy $M/README_tstar1.30.txt claw:$C/
ssh claw "sudo chown -R www-data:www-data $C"
for f in ${P}_wabs.npy ${P}_wz.npy ${P}_lam2.npy; do
  [ "$(ssh claw "sudo stat -c %s $C/$f")" = "$(stat -f %z $M/$f)" ] && [ "$(stat -f %z $M/$f)" = 268435584 ] || { echo "size check failed $f"; exit 1; }
done
echo "all on the Mac and claw"; ssh claw "ls -la $C"
rm -f $M/l40kbase_m512.slab; ssh hal "rm -rf $H/r40kf $H/l40kbase_m512.slab"; ssh rc "rm -f $Y/l40kbase_m512.slab"
