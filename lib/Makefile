CUDA	:= $(CUDA_INSTALL_PATH)
INC	:= -I$(CUDA)/include -I$(CUDA)/samples/common/inc/ -I.
LIB	:= -L$(CUDA)/lib64
BUILD_DIR	:=../build
ASTROLIB_DIR :=../astrolib

# CUDA code generation flags
GENCODE_SM20    := -gencode arch=compute_20,code=sm_20 # Fermi
GENCODE_SM21    := -gencode arch=compute_21,code=sm_21 # Fermi
GENCODE_SM30    := -gencode arch=compute_30,code=sm_30 # Kepler
GENCODE_SM32    := -gencode arch=compute_32,code=sm_32 # Kepler
GENCODE_SM35	:= -gencode arch=compute_35,code=sm_35 # Kepler
GENCODE_SM37	:= -gencode arch=compute_37,code=sm_37 # Kepler
GENCODE_SM50	:= -gencode arch=compute_50,code=sm_50 # Maxwell
GENCODE_SM52	:= -gencode arch=compute_52,code=sm_52 # Maxwell
GENCODE_SM53	:= -gencode arch=compute_53,code=sm_53 # Maxwell
GENCODE_SM60	:= -gencode arch=compute_60,code=sm_60 # Pascal
GENCODE_SM61	:= -gencode arch=compute_61,code=sm_61 # Pascal
GENCODE_FLAGS   := $(GENCODE_SM60)

ifeq ($(cache),off)
        NVCCFLAGS := $(INC) -g -lineinfo -Xcompiler -fopenmp -O3 -lm --use_fast_math\
         --ptxas-options=-v -Xptxas -dlcm=cg $(GENCODE_FLAGS)
else
        NVCCFLAGS := $(INC) -g -lineinfo -Xcompiler -fopenmp -O3 -lm --use_fast_math\
         --ptxas-options=-v -lcuda -lcudart  -lcurand -lcufft -lcudadevrt -Xptxas -dlcm=cg $(GENCODE_FLAGS)
endif

#ifeq ($(gpu),sm_35)
#        NVCCFLAGS += -D SM_35
#endif


all:	${BUILD_DIR} $(BUILD_DIR)/dedisperse-gpu ${ASTROLIB_DIR} ${ASTROLIB_DIR}/libastrolib.a

$(BUILD_DIR)/%.o : %.cu
	nvcc -o $@ -c $(NVCCFLAGS) $<

$(BUILD_DIR)/%.o : %.cpp
	nvcc -o $@ -c $(NVCCFLAGS) $<
	
$(BUILD_DIR)/%.o : %.c
	nvcc -o $@ -c $(NVCCFLAGS) $<
	
${ASTROLIB_DIR}:
	mkdir -p ${ASTROLIB_DIR}

${BUILD_DIR}:
	mkdir -p ${BUILD_DIR}	

$(BUILD_DIR)/main.o : 					headers/headers_mains.h
$(BUILD_DIR)/host_main_function.o:		headers/host_main_function.h 
$(BUILD_DIR)/host_allocate_memory.o:	headers/host_allocate_memory.h
$(BUILD_DIR)/host_acceleration.o:		headers/host_acceleration.h
$(BUILD_DIR)/host_analysis.o:			headers/host_analysis.h 
$(BUILD_DIR)/host_periods.o:			headers/host_periods.h 
$(BUILD_DIR)/host_debug.o:	 			headers/host_debug.h
$(BUILD_DIR)/host_get_file_data.o: 		headers/host_get_file_data.h
$(BUILD_DIR)/host_get_user_input.o:	 	headers/host_get_user_input.h
$(BUILD_DIR)/host_get_recorded_data.o:	headers/host_get_recorded_data.h
$(BUILD_DIR)/host_help.o:	 			headers/host_help.h
$(BUILD_DIR)/host_rfi.o:	 			headers/host_rfi.h
$(BUILD_DIR)/host_stratagy.o:	 		headers/host_stratagy.h
$(BUILD_DIR)/host_statistics.o:	 		headers/host_statistics.h
$(BUILD_DIR)/host_write_file.o:	 		headers/host_write_file.h


#$(BUILD_DIR)/device_main.o:	device_main.cu
#	nvcc -c $(NVCCFLAGS) device_main.cu

$(ASTROLIB_DIR)/libastrolib.a: $(BUILD_DIR)/host_acceleration.o $(BUILD_DIR)/host_write_file.o $(BUILD_DIR)/device_main.o $(BUILD_DIR)/host_allocate_memory.o \
	$(BUILD_DIR)/host_analysis.o $(BUILD_DIR)/host_periods.o $(BUILD_DIR)/host_debug.o $(BUILD_DIR)/host_get_file_data.o $(BUILD_DIR)/host_get_user_input.o \
	$(BUILD_DIR)/host_get_recorded_data.o $(BUILD_DIR)/host_help.o $(BUILD_DIR)/host_rfi.o $(BUILD_DIR)/host_stratagy.o $(BUILD_DIR)/host_statistics.o \
	$(BUILD_DIR)/host_main_function.o $(BUILD_DIR)/fresnl.o $(BUILD_DIR)/median.o $(BUILD_DIR)/presto_funcs.o $(BUILD_DIR)/characteristics.o $(BUILD_DIR)/dcdflib.o $(BUILD_DIR)/ipmpar.o \
	$(BUILD_DIR)/fdas_host.o  $(BUILD_DIR)/fdas_device.o $(BUILD_DIR)/fdas_util.o
	ar rcs $(ASTROLIB_DIR)/libastrolib.a $(BUILD_DIR)/host_acceleration.o $(BUILD_DIR)/host_write_file.o $(BUILD_DIR)/device_main.o $(BUILD_DIR)/host_allocate_memory.o \
	$(BUILD_DIR)/host_analysis.o $(BUILD_DIR)/host_periods.o $(BUILD_DIR)/host_debug.o $(BUILD_DIR)/host_get_file_data.o $(BUILD_DIR)/host_get_user_input.o \
	$(BUILD_DIR)/host_get_recorded_data.o $(BUILD_DIR)/host_help.o $(BUILD_DIR)/host_rfi.o $(BUILD_DIR)/host_stratagy.o $(BUILD_DIR)/host_statistics.o \
	$(BUILD_DIR)/host_main_function.o $(BUILD_DIR)/fresnl.o $(BUILD_DIR)/median.o $(BUILD_DIR)/presto_funcs.o $(BUILD_DIR)/characteristics.o $(BUILD_DIR)/dcdflib.o $(BUILD_DIR)/ipmpar.o \
	$(BUILD_DIR)/fdas_host.o  $(BUILD_DIR)/fdas_device.o $(BUILD_DIR)/fdas_util.o

$(BUILD_DIR)/dedisperse-gpu: $(BUILD_DIR)/main.o $(ASTROLIB_DIR)/libastrolib.a
	nvcc -o $(BUILD_DIR)/dedisperse-gpu $(BUILD_DIR)/main.o -L$(ASTROLIB_DIR)/ -lastrolib -L${LIB} $(NVCCFLAGS)

clean:
	rm -f $(BUILD_DIR)/dedisperse-gpu $(BUILD_DIR)/*.a $(BUILD_DIR)/*.o $(ASTROLIB_DIR)/*.a
