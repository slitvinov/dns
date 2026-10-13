NVCC = nvcc
NVCCFLAGS = -O3 -arch=sm_90
LDLIBS = -lcufft

all: fourier tg etd box boxs boxsym symq symfull sp2 vm2 vm2_g8.bin

fourier: fourier.cu
	$(NVCC) $(NVCCFLAGS) -o fourier fourier.cu $(LDLIBS)

tg: tg.cu
	$(NVCC) $(NVCCFLAGS) -o tg tg.cu $(LDLIBS)

etd: etd.cu
	$(NVCC) $(NVCCFLAGS) -o etd etd.cu $(LDLIBS)

box: box.cu
	$(NVCC) $(NVCCFLAGS) -o box box.cu $(LDLIBS)

boxs: boxs.cu
	$(NVCC) $(NVCCFLAGS) -o boxs boxs.cu $(LDLIBS)

boxsym: boxsym.cu
	$(NVCC) $(NVCCFLAGS) -o boxsym boxsym.cu $(LDLIBS)

vm3: vm3.cu
	$(NVCC) $(NVCCFLAGS) -o vm3 vm3.cu $(LDLIBS)

sp2: sp2.cu
	$(NVCC) $(NVCCFLAGS) -o sp2 sp2.cu $(LDLIBS)

vm2: vm2.cu
	$(NVCC) $(NVCCFLAGS) -default-stream per-thread -Xcompiler -ffp-contract=off -o vm2 vm2.cu $(LDLIBS)

vm2_g8.bin: gen_vm2.py
	python3 gen_vm2.py

symq: symq.c
	$(CC) -O2 -o symq symq.c

symfull: symfull.c
	$(CC) -O2 -o symfull symfull.c

clean:
	rm -f fourier tg etd box boxs boxsym symq symfull sp2 vm2 vm2_g4.bin vm2_g8.bin vm3
