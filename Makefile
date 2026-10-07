CFLAGS = -O3 -march=native -fopenmp
LDLIBS = -lfftw3 -lfftw3_omp -lm

all: fourier tg

fourier: fourier.c
	$(CC) $(CFLAGS) -o fourier fourier.c $(LDLIBS)

tg: tg.c
	$(CC) $(CFLAGS) -o tg tg.c $(LDLIBS)

clean:
	rm -f fourier tg
