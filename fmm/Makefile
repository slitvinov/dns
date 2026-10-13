NVCC = nvcc
NVCCFLAGS = -O3 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_90,code=sm_90

all: fmm2 gen2 fmm3 gen3

fmm2: fmm2.cu
	$(NVCC) $(NVCCFLAGS) -o fmm2 fmm2.cu

gen2: gen2.c
	$(CC) -O2 -o gen2 gen2.c -lm

fmm3: fmm3.cu
	$(NVCC) $(NVCCFLAGS) -o fmm3 fmm3.cu -lcublas -lcufft

gen3: gen3.c
	$(CC) -O2 -o gen3 gen3.c -lm
