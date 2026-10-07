NVCC = nvcc
NVCCFLAGS = -O3 -arch=sm_90
LDLIBS = -lcufft

all: fourier tg etd

fourier: fourier.cu
	$(NVCC) $(NVCCFLAGS) -o fourier fourier.cu $(LDLIBS)

tg: tg.cu
	$(NVCC) $(NVCCFLAGS) -o tg tg.cu $(LDLIBS)

etd: etd.cu
	$(NVCC) $(NVCCFLAGS) -o etd etd.cu $(LDLIBS)

clean:
	rm -f fourier tg etd
