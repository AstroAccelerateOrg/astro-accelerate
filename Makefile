CUDA	:= $(CUDA_INSTALL_PATH)
CUDASDK := -I/usr/local/cudaSDK/C/common/inc/ -I.
INC	:= -I$(CUDA)/include -I.
LIB	:= -L$(CUDA)/lib64
LIBS	:= -L$(CUDA)/lib

ifeq ($(cache),off)
        NVCCFLAGS := -O3 -arch=sm_21 --ptxas-options=-v --use_fast_math -Xptxas -dlcm=cg
else
        NVCCFLAGS := -O3 -arch=sm_21 --ptxas-options=-v --use_fast_math
endif

ICCFLAGS := -i-static -fPIC -xSSE4.1 -O3 -parallel -openmp -par-report -vec-report

all:	dedisperse-gpu

host_main.o:	host_main.c headers_mains.h 
	icc host_main.c -c -o host_main.o $(ICCFLAGS) $(INC) $(LIB) $(CUDASDK)

host_analysis.o:	host_analysis.c host_analysis.h 
	icc host_analysis.c -c -o host_analysis.o $(ICCFLAGS) $(INC) $(LIB) $(CUDASDK)

host_bin.o:	host_bin.c host_bin.h 
	icc host_bin.c -c -o host_bin.o $(ICCFLAGS) $(INC) $(LIB) $(CUDASDK)

host_init.o:	host_init.c host_init.h
	icc host_init.c -c -o host_init.o $(ICCFLAGS) $(INC) $(LIB) $(CUDASDK)

host_help.o:	host_help.c host_help.h
	icc host_help.c -c -o host_help.o $(ICCFLAGS) $(INC) $(LIB) $(CUDASDK)

device_main.o:	device_main.cu
	nvcc device_main.cu -c -o device_main.o $(INC) $(CUDASDK) $(LIB) $(NVCCFLAGS)

dedisperse-gpu:	host_main.o host_analysis.o host_bin.o host_init.o host_help.o device_main.o 
	icc -fPIC -o dedisperse-gpu host_main.o host_bin.o host_analysis.o host_init.o host_help.o device_main.o $(ICCFLAGS) $(LIB) $(INC) $(LIB) $(CUDASDK) $(LIBS) -lcudart

clean:
	rm -f dedisperse-gpu *.a *.o
