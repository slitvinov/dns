<h2>Build</h2>

<pre>
$ mpicc mpi/main.c -lfftw3_mpi -lfftw3 -O2 -g -lm
$ mpiexec ./a.out
k = 1.2499062522321387e-01
k = 1.2498125036872199e-01
k = 1.2497187465530142e-01
k = 1.2496249730189801e-01
k = 1.2495311752707980e-01
</pre>
