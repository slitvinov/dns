NVCC = nvcc
NVCCFLAGS = -O3 -arch=sm_90
LDLIBS = -lcufft

all: fourier tg etd box

fourier: fourier.cu
	$(NVCC) $(NVCCFLAGS) -o fourier fourier.cu $(LDLIBS)

tg: tg.cu
	$(NVCC) $(NVCCFLAGS) -o tg tg.cu $(LDLIBS)

etd: etd.cu
	$(NVCC) $(NVCCFLAGS) -o etd etd.cu $(LDLIBS)

box: box.cu
	$(NVCC) $(NVCCFLAGS) -o box box.cu $(LDLIBS)

clean:
	rm -f fourier tg etd box
