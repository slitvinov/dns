CC = mpicc
CFLAGS = -O3 -march=native -fopenmp
LDLIBS = -lfftw3_mpi -lfftw3_omp -lfftw3 -lm

all: fourier tg

fourier: fourier.c
	$(CC) $(CFLAGS) -o fourier fourier.c $(LDLIBS)

tg: tg.c
	$(CC) $(CFLAGS) -o tg tg.c $(LDLIBS)

clean:
	rm -f fourier tg
