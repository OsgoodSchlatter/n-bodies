CC	= gcc
CFLAGS	= -O0 -g -Wall
LDFLAGS = -lm -lpthread -lX11
VERBOSE	=
TARGET	= nbody_brute_force #nbody_barnes_hut
OBJS	= ui.o xstuff.o nbody_tools.o nbody_alloc.o

CUDA_HOME=/usr/local/cuda
CUDASDK=/usr/local/cuda/samples
CUDANVIDIA=/usr/lib64/nvidia
NVCC=$(CUDA_HOME)/bin/nvcc -arch=sm_35

NVCFLAGS= -I $(CUDA_HOME)/include -I$(CUDASDK)/common/inc
CFLAGS= -L $(CUDA_HOME)/lib64 -L $(CUDANVIDIA)

DISPLAY = -DDISPLAY
#DISPLAY =

#DUMP = -DDUMP_RESULT
DUMP =

all: $(TARGET)

cuda_recurence: test_recurence.o
	$(NVCC) -arch=sm_35 $(CFLAGS)  -o $@ $^

#nbody_brute_force: nbody_brute_force.o
#	$(NVCC) -arch=sm_35 $(CFLAGS)  -o $@ $^

#nbody_barnes_hut: nbody_barnes_hut.cu $(OBJS)
#	$(NVCC) $(VERBOSE) -o $@ $< $(OBJS)  $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< $(VERBOSE) $(DISPLAY) $(DUMP)

%.o: %.cu
	$(NVCC) -arch=sm_35 $(NVCFLAGS) -c -o $@ $<

clean:
	rm -f *.o $(TARGET)
