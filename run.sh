set 100 200 400 800 1600 3000

for i
do echo $i
done | awk '{printf "%04d %.16e\n", $1, 1/$1}' | \
    xargs --process-slot-var I -n 2 -P `nproc` sh -c \
	  'echo taskset --cpu-list $i -- -T 10 -n $1 -i tgv.raw $0' 
