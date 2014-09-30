CUDA	:= $(CUDA_INSTALL_PATH)
INC	:= -I$(CUDA)/include -I$(CUDA)/samples/common/inc/ -I.
LIB	:= -L$(CUDA)/lib64

ifeq ($(cache),off)
        NVCCFLAGS := -g -Xcompiler -fopenmp -O3 -lm -arch=sm_35 --use_fast_math --ptxas-options=-v -Xptxas -dlcm=cg
else
        NVCCFLAGS := -g -Xcompiler -fopenmp -O3 -lm -arch=sm_35 --use_fast_math --ptxas-options=-v -lcuda -lcudart  -lcurand -lcufft
endif

all:	dedisperse-gpu

host_main.o:	host_main.cu headers_mains.h 
	nvcc -c $(NVCCFLAGS) host_main.cu 

host_allocate_memory.o:	host_allocate_memory.cu host_allocate_memory.h
	nvcc -c $(NVCCFLAGS) host_allocate_memory.cu 

host_analysis.o:	host_analysis.cu host_analysis.h 
	nvcc -c $(NVCCFLAGS) host_analysis.cu

host_periods.o:	host_periods.cu host_periods.h 
	nvcc -c $(NVCCFLAGS) host_periods.cu

host_debug.o:	host_debug.cu host_debug.h
	nvcc -c $(NVCCFLAGS) host_debug.cu 

host_get_file_data.o:	host_get_file_data.cu host_get_file_data.h
	nvcc -c $(NVCCFLAGS) host_get_file_data.cu

host_get_user_input.o:	host_get_user_input.cu host_get_user_input.h
	nvcc -c $(NVCCFLAGS) host_get_user_input.cu

host_get_recorded_data.o:	host_get_recorded_data.cu host_get_recorded_data.h
	nvcc -c $(NVCCFLAGS) host_get_recorded_data.cu 

host_help.o:	host_help.cu host_help.h
	nvcc -c $(NVCCFLAGS) host_help.cu

host_rfi.o:	host_rfi.cu host_rfi.h
	nvcc -c $(NVCCFLAGS) host_rfi.cu

host_stratagy.o:	host_stratagy.cu host_stratagy.h
	nvcc -c $(NVCCFLAGS) host_stratagy.cu

host_statistics.o:	host_statistics.cu host_statistics.h
	nvcc -c $(NVCCFLAGS) host_statistics.cu

host_write_file.o:	host_write_file.cu host_write_file.h
	nvcc -c $(NVCCFLAGS) host_write_file.cu

device_main.o:	device_main.cu
	nvcc -c $(NVCCFLAGS) device_main.cu

dedisperse-gpu:	host_main.o host_allocate_memory.o host_analysis.o host_periods.o host_debug.o host_get_file_data.o host_get_user_input.o host_get_recorded_data.o host_help.o host_rfi.o host_stratagy.o host_statistics.o host_write_file.o device_main.o
	nvcc -o dedisperse-gpu $(NVCCFLAGS) host_main.o host_allocate_memory.o host_analysis.o host_periods.o host_debug.o host_get_file_data.o host_get_user_input.o host_get_recorded_data.o host_help.o host_rfi.o host_stratagy.o host_statistics.o host_write_file.o device_main.o

clean:
	rm -f dedisperse-gpu *.a *.o
