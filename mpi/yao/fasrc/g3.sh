#!/bin/bash
# g3.sh: Re 9000 gates for the Re 40000 run, then the Re 9000 video fields (4 H200, one node)
#  R1: 512 to t* 1.3, k. dump, padded to 1024 to t 5.6 (compare with r9k1024);  R2: 768 without -S (compare with r9k768)
. /etc/profile > /dev/null 2>&1
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; . $Y/src/mpenv.sh; B=$Y/src/boxsymmp_b4
PY=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/venv/bin/python
mp() { mpirun -np 4 --oversubscribe --bind-to none "$@"; }
G="-X 2 -Y 2 -Z 2 -n 8.533761e-4 -s 0.01 -f -c 1"
mkdir -p $Y/run/r9kpad && cd $Y/run/r9kpad
mp $B -x 512 -y 512 -z 512 $G -S -k -I yao -T 2.7911 -t 2.7911 > out.0 2> err.0
k=$(ls k.[0-9]*[0-9] | tail -1)
mp $B -x 1024 -y 1024 -z 1024 $G -S -e 0.25 -r $k -G 512,512,512 -t 5.6 > out.1 2> err.1; rm -f k.*; touch done
mkdir -p $Y/run/r9k768noS && cd $Y/run/r9k768noS
mp $B -x 768 -y 768 -z 768 $G -e 0.25 -I yao -t 5.6 > out 2> err; touch done
mkdir -p $Y/run/r9kvid/fields && cd $Y/run/r9kvid
mp $B -x 1024 -y 1024 -z 1024 $G -S -k -e 0.25 -I yao -d 100 -T 0.5367436332,1.073487266,1.6102309,2.146974533,2.683718166,3.220461799,3.757205432,4.293949065,4.830692699,5.367436332 -t 5.367436332 > out 2> err
for k in k.[0-9]*[0-9]; do ts=$(awk '{printf "%.2f", $1 * 0.4657717102700246}' $k.t); OMP_NUM_THREADS=16 $PY $Y/src/k2fields.py $k 1024 512 fields/re9000_tstar$ts >> conv.log 2>&1 && rm -f $k; done
touch done
# Re 5000 at 1024 (C 1, -S, to t* 2.6): the |w|max convergence 512 -> 768 -> 1024 against vm3 2pi/2048 window
mkdir -p $Y/run/r5k1024 && cd $Y/run/r5k1024
mp $B -x 1024 -y 1024 -z 1024 -X 2 -Y 2 -Z 2 -n 1.5360770488e-3 -s 0.01 -f -c 1 -S -e 0.25 -I yao -t 5.6 > out 2> err; touch done
scancel $SLURM_JOB_ID
