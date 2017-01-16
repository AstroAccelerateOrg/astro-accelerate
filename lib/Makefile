CUDA	:= $(CUDA_INSTALL_PATH)
INC	:= -I$(CUDA)/include -I$(CUDA)/samples/common/inc/ -I.
LIB	:= -L$(CUDA)/lib64
BUILD_DIR	:=../build
ASTROLIB_DIR :=../astrolib

CCFLAGS = $(INC) -g -lineinfo -O3 -std=c++11

CCFLAGS += -DSM_35

ifeq ($(cache),off)
        NVCCFLAGS := $(CCFLAGS) -Xcompiler -fopenmp -arch=compute_35 --use_fast_math\
         --ptxas-options=-v -Xptxas -dlcm=cg
else
        NVCCFLAGS := $(CCFLAGS) -Xcompiler -fopenmp -arch=compute_35 --use_fast_math\
         --ptxas-options=-v -lcuda -lcudart -lnppc -lnppi -lcurand -lcufft -Xptxas -dlcm=cg
endif

NVCCFLAGS += -gencode arch=compute_35,code=sm_35
NVCCFLAGS += -gencode arch=compute_37,code=sm_37
NVCCFLAGS += -gencode arch=compute_52,code=sm_52
#NVCCFLAGS += -gencode arch=compute_60,code=sm_60
#NVCCFLAGS += -gencode arch=compute_61,code=sm_61


NVCCFLAGS += --source-in-ptx --keep


all: $(BUILD_DIR)/dedisperse-gpu ${BUILD_DIR}/peak_find_test ${ASTROLIB_DIR}/libastrolib.a

$(BUILD_DIR)/host_%.o: host_%.cu
	$(CXX) -ansi -fopenmp -x c++ -o $@ -c -Wall $(CCFLAGS) $<

$(BUILD_DIR)/%.o : %.cu
	$(CUDA)/bin/nvcc -o $@ -c $(NVCCFLAGS) $<

$(BUILD_DIR)/%.o : %.cpp
	$(CUDA)/bin/nvcc -o $@ -c $(NVCCFLAGS) $<

${ASTROLIB_DIR}:
	mkdir -p ${ASTROLIB_DIR}

${BUILD_DIR}:
	mkdir -p ${BUILD_DIR}


$(BUILD_DIR)/main.o : 					AstroAccelerate/headers_mains.h

$(BUILD_DIR)/host_main_function.o:		AstroAccelerate/host_main_function.h

$(BUILD_DIR)/host_allocate_memory.o:	AstroAccelerate/host_allocate_memory.h

$(BUILD_DIR)/host_acceleration.o:		AstroAccelerate/host_acceleration.h

$(BUILD_DIR)/host_analysis.o:			AstroAccelerate/host_analysis.h

$(BUILD_DIR)/host_periods.o:			AstroAccelerate/host_periods.h

$(BUILD_DIR)/host_debug.o:	 			AstroAccelerate/host_debug.h

$(BUILD_DIR)/host_get_file_data.o: 		AstroAccelerate/host_get_file_data.h

$(BUILD_DIR)/host_get_user_input.o:	 	AstroAccelerate/host_get_user_input.h

$(BUILD_DIR)/host_get_recorded_data.o:	AstroAccelerate/host_get_recorded_data.h

$(BUILD_DIR)/host_help.o:	 			AstroAccelerate/host_help.h

$(BUILD_DIR)/host_rfi.o:	 			AstroAccelerate/host_rfi.h

$(BUILD_DIR)/host_stratagy.o:	 		AstroAccelerate/host_stratagy.h

$(BUILD_DIR)/host_statistics.o:	 		AstroAccelerate/host_statistics.h

$(BUILD_DIR)/host_write_file.o:	 		AstroAccelerate/host_write_file.h

#$(BUILD_DIR)/device_main.o:	device_main.cu
#	nvcc -c $(NVCCFLAGS) device_main.cu

$(ASTROLIB_DIR)/libastrolib.a: $(BUILD_DIR)/host_acceleration.o $(BUILD_DIR)/host_write_file.o $(BUILD_DIR)/device_main.o $(BUILD_DIR)/host_allocate_memory.o \
	$(BUILD_DIR)/host_analysis.o $(BUILD_DIR)/host_periods.o $(BUILD_DIR)/host_debug.o $(BUILD_DIR)/host_get_file_data.o $(BUILD_DIR)/host_get_user_input.o \
	$(BUILD_DIR)/host_get_recorded_data.o $(BUILD_DIR)/host_help.o $(BUILD_DIR)/host_rfi.o $(BUILD_DIR)/host_stratagy.o $(BUILD_DIR)/host_statistics.o \
	$(BUILD_DIR)/host_main_function.o
	ar rcs $(ASTROLIB_DIR)/libastrolib.a $(BUILD_DIR)/host_acceleration.o $(BUILD_DIR)/host_write_file.o $(BUILD_DIR)/device_main.o $(BUILD_DIR)/host_allocate_memory.o \
	$(BUILD_DIR)/host_analysis.o $(BUILD_DIR)/host_periods.o $(BUILD_DIR)/host_debug.o $(BUILD_DIR)/host_get_file_data.o $(BUILD_DIR)/host_get_user_input.o \
	$(BUILD_DIR)/host_get_recorded_data.o $(BUILD_DIR)/host_help.o $(BUILD_DIR)/host_rfi.o $(BUILD_DIR)/host_stratagy.o $(BUILD_DIR)/host_statistics.o \
	$(BUILD_DIR)/host_main_function.o

$(BUILD_DIR)/dedisperse-gpu: $(BUILD_DIR)/main.o $(ASTROLIB_DIR)/libastrolib.a
	$(CUDA)/bin/nvcc -o $(BUILD_DIR)/dedisperse-gpu $(BUILD_DIR)/main.o -L$(ASTROLIB_DIR)/ -lastrolib ${LIB} $(NVCCFLAGS)

$(BUILD_DIR)/peak_find_test: $(BUILD_DIR)/host_peak_find_main.o $(ASTROLIB_DIR)/libastrolib.a
	$(CUDA)/bin/nvcc -o $(BUILD_DIR)/peak_find_test $(BUILD_DIR)/host_peak_find_main.o -L$(ASTROLIB_DIR)/ -lastrolib ${LIB} $(NVCCFLAGS)

clean:
	rm -f $(BUILD_DIR)/dedisperse-gpu $(BUILD_DIR)/*.a $(BUILD_DIR)/*.o $(ASTROLIB_DIR)/*.a
