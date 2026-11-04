NVCC = nvcc
NVCCFLAGS = -O3 -arch=sm_90
LDLIBS = -lcufft

all: fourier tg tgr rview etd box boxs boxsym symq symfull

fourier: fourier.cu
	$(NVCC) $(NVCCFLAGS) -o fourier fourier.cu $(LDLIBS)

tg: tg.cu
	$(NVCC) $(NVCCFLAGS) -o tg tg.cu $(LDLIBS)

tgr: tgr.cu render.h
	$(NVCC) $(NVCCFLAGS) -o tgr tgr.cu $(LDLIBS)

rview: rview.cu render.h
	$(NVCC) $(NVCCFLAGS) -o rview rview.cu

etd: etd.cu
	$(NVCC) $(NVCCFLAGS) -o etd etd.cu $(LDLIBS)

box: box.cu
	$(NVCC) $(NVCCFLAGS) -o box box.cu $(LDLIBS)

boxs: boxs.cu
	$(NVCC) $(NVCCFLAGS) -o boxs boxs.cu $(LDLIBS)

boxsym: boxsym.cu
	$(NVCC) $(NVCCFLAGS) -o boxsym boxsym.cu $(LDLIBS)

symq: symq.c
	$(CC) -O2 -o symq symq.c

symfull: symfull.c
	$(CC) -O2 -o symfull symfull.c

clean:
	rm -f fourier tg tgr rview etd box boxs boxsym symq symfull
