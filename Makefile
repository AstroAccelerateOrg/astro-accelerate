################################################################################
# Makefile for astro-accelerate
#
# Description: Custom Makefile that can be used instead of CMakeLists.txt
#
################################################################################
CUDA	:= $(CUDA_INSTALL_PATH)
INC	:= -I$(CUDA)include -I$(CUDA)samples/common/inc/
LIB	:= -L$(CUDA)/lib64
BUILD_DIR	:=./obj/
ASTROLIB_DIR	:=./
SRC_DIR		:= src/

CPP_FILES := $(wildcard src/*.cpp)
CU_FILES := ${wildcard src/*.cu}
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o))) $(addprefix obj/,$(notdir $(CU_FILES:.cu=.o)))
EXAMPLES_FILES := examples_dedispersion examples_periodicity examples_dedispersion_and_analysis examples_filterbank_dedispersion

CXXFLAGS	:= -O3
LDFLAGS		= `root-config --libs`
COMPILEJOBS	= astro-accelerate
HEADERS		= include/*.hpp include/*.cuh
INCLUDE		= -Iinclude

# CUDA code generation flags
GENCODE_SM35	:= -gencode arch=compute_35,code=sm_35 # Kepler
GENCODE_SM37	:= -gencode arch=compute_37,code=sm_37 # Kepler
GENCODE_SM50	:= -gencode arch=compute_50,code=sm_50 # Maxwell
GENCODE_SM52	:= -gencode arch=compute_52,code=sm_52 # Maxwell
GENCODE_SM53	:= -gencode arch=compute_53,code=sm_53 # Maxwell
GENCODE_SM60	:= -gencode arch=compute_60,code=sm_60 # Pascal
GENCODE_SM61	:= -gencode arch=compute_61,code=sm_61 # Pascal
GENCODE_SM70	:= -gencode arch=compute_70,code=sm_70 # Volta
GENCODE_SM75    := -gencode arch=compute_75,code=sm_75 # Turing
GENCODE_FLAGS   := $(GENCODE_SM61)

ifeq ($(cache),off)
        NVCCFLAGS := $(INC) ${INCLUDE} -g -lineinfo -Xcompiler -O3 -lm --use_fast_math\
        --ptxas-options=-v -Xptxas -dlcm=cg $(GENCODE_FLAGS)
else
        NVCCFLAGS := $(INC) ${INCLUDE} -g -lineinfo -Xcompiler -O3 -lm --use_fast_math\
        --ptxas-options=-v -lcuda -lcudart  -lcurand -lcufft -lcudadevrt -Xptxas -dlcm=cg $(GENCODE_FLAGS)
endif

LIBJOBS := libastroaccelerate.a

all:	MAKE_OBJ_FOLDER ${LIBJOBS} ${COMPILEJOBS} ${EXAMPLES_FILES}

$(BUILD_DIR)%.o :	${SRC_DIR}%.cu
			@echo Compiling $@ ...
			nvcc $(NVCCFLAGS) -c -o $@ $<

$(BUILD_DIR)%.o :	${SRC_DIR}%.cpp
			@echo Compiling $@ ...
			nvcc $(NVCCFLAGS) -c -o $@ $<

MAKE_OBJ_FOLDER :
			mkdir -p obj/

libastroaccelerate.a: ${OBJ_FILES}
			@echo Making static library libastroaccelerate.a
			ar rcs libastroaccelerate.a $(OBJ_FILES)

astro-accelerate: libastroaccelerate.a
			nvcc -o astro-accelerate $(OBJ_FILES) -L$(ASTROLIB_DIR) -lastroaccelerate -L${LIB} $(NVCCFLAGS)

examples_dedispersion: libastroaccelerate.a
			nvcc -o examples_dedispersion ./examples/src/dedispersion.cpp -L$(ASTROLIB_DIR) -lastroaccelerate -L${LIB} $(NVCCFLAGS)

examples_periodicity: libastroaccelerate.a
			nvcc -o examples_periodicity ./examples/src/periodicity.cpp -L$(ASTROLIB_DIR) -lastroaccelerate -L${LIB} $(NVCCFLAGS)

examples_dedispersion_and_analysis: libastroaccelerate.a
			nvcc -o examples_dedispersion_and_analysis ./examples/src/dedispersion_and_analysis.cpp -L$(ASTROLIB_DIR) -lastroaccelerate -L${LIB} $(NVCCFLAGS)

examples_filterbank_dedispersion: libastroaccelerate.a
			nvcc -o examples_filterbank_dedispersion ./examples/src/filterbank_dedispersion.cpp -L$(ASTROLIB_DIR) -lastroaccelerate -L${LIB} $(NVCCFLAGS)

clean:
	rm -f astro-accelerate *.a $(BUILD_DIR)*.o $(ASTROLIB_DIR)*.a
