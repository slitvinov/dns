<h2>Build</h2>

Serial
<pre>
$ c99 main.c -lfftw3_mpi -lfftw3 -O2 -g -lm
$ ./a.out
eng = 1.2500000000000000e-01
eng = 1.2499062522136267e-01
eng = 1.2498125036475010e-01
eng = 1.2497187464920783e-01
eng = 1.2496249729368313e-01
eng = 1.2495311751674273e-01
</pre>

with MPI
<pre>
$ mpicc mpi/main.c -lfftw3_mpi -lfftw3 -O2 -g -lm
$ mpiexec ./a.out
eng: 1.2499062522321387e-01
eng: 1.2498125036872199e-01
eng: 1.2497187465530142e-01
eng: 1.2496249730189801e-01
eng: 1.2495311752707980e-01
</pre>
