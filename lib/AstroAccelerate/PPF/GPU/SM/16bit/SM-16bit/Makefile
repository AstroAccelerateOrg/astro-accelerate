###############################################################
# Makefile for Artemis project -- Polyphase filter
# created by JaN
# contact: jan.novotny@fpf.slu.cz, karel.adamek@fpf.slu.cz,
#	   wes.armour@oerc.ox.ac.uk
###############################################################

###############################################################
# CUDA_HOME are supposed to be on default position
# and set it in your PATH .bashrc
###############################################################
SDK := /home/novotny/NVIDIA_GPU_Computing_SDK/CUDALibraries/common/inc/
INC := -I/usr/local/cuda/include #-I$(SDK)
LIB := -L/usr/local/cuda/lib64 -lcudart -lcufft -lfftw3f -lcuda

# use this compilers
# g++ just because the file write
GCC = g++
NVCC = nvcc


###############################################################
# Basic flags for compilers, one for debug options
# fmad flags used for reason of floating point operation
###############################################################
NVCCFLAGS = -O2 -arch=sm_52 --ptxas-options=-v --use_fast_math -Xcompiler -Wextra -lineinfo
#NVCCFLAGS = -O0 -arch=sm_20 --ptxas-options=-v -Xptxas -dlcm=cg --use_fast_math -Xcompiler -Wextra
#NVCCFLAGS= -g -G -arch=sm_20 --ptxas-options=-v --use_fast_math -Xcompiler -Wextra
GCC_OPTS =-O3 -Wall -Wextra $(INC)
ANALYZE = SM-polyphase.exe

all: clean analyze

analyze: polyphase.o SM-maxwell-16bit.o reference.o Makefile
	$(NVCC) -o $(ANALYZE) SM-maxwell-16bit.o reference.o polyphase.o $(LIB) $(NVCCFLAGS) 

SM-maxwell-16bit.o: ../timer.h ../utils_cuda.h ../utils_file.h
	$(NVCC) -c SM-maxwell-16bit.cu $(NVCCFLAGS)

polyphase.o: ../polyphase.c
	$(GCC) -c ../polyphase.c $(GCC_OPTS)

reference.o: ../reference.c 
	$(GCC) -c ../reference.c $(GCC_OPTS) $(INC) -lfftw3

clean:	
	rm -f *.o *.~ $(ANALYZE)


