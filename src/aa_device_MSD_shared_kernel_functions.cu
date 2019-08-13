#include "aa_device_MSD_shared_kernel_functions.cuh"
#include "aa_device_MSD_shared_kernel_functions.hpp"

namespace astroaccelerate {

  // Computes mean and standard deviation from partial
  __global__ void MSD_GPU_final_regular(float *d_input, float *d_output, int size) {
    __shared__ float s_input[3*WARP*WARP];

    float M, S, j;

    Sum_partials_regular( &M, &S, &j, d_input, s_input, size);

    //----------------------------------------------
    //---- Writing data
    if (threadIdx.x == 0) {
      d_output[0] = M / j;
      d_output[1] = sqrt(S / j);
      d_output[2] = j;
    }
  }

  __global__ void MSD_GPU_final_regular(float *d_input, float *d_MSD, float *d_pp, int size) {
    __shared__ float s_input[3*WARP*WARP];

    float M, S, j;

    Sum_partials_regular( &M, &S, &j, d_input, s_input, size);

    if(d_pp[2]>0){
      //Merge(&M, &S, &j, d_pp[0]*d_pp[2], (d_pp[1]*d_pp[1])*d_pp[2], d_pp[2]);
      Merge(&M, &S, &j, d_pp[0], d_pp[1], d_pp[2]);
    }

    //----------------------------------------------
    //---- Writing data
    if (threadIdx.x == 0) {
      d_MSD[0] = M / j;
      d_MSD[1] = sqrt(S / j);
      d_MSD[2] = j;
      d_pp[0] = M;
      d_pp[1] = S;
      d_pp[2] = j;
    }
  }

  __global__ void MSD_GPU_final_nonregular(float *d_input, float *d_MSD, int size) {
    __shared__ float s_input[3*WARP*WARP];

    float M, S, j;

    Sum_partials_nonregular( &M, &S, &j, d_input, s_input, size);

    //----------------------------------------------
    //---- Writing data
    if (threadIdx.x == 0) {
      d_MSD[0] = M / j;
      d_MSD[1] = sqrt(S / j);
      d_MSD[2] = j;
    }
  }

  __global__ void MSD_GPU_final_nonregular(float *d_input, float *d_MSD, float *d_pp, int size) {
    __shared__ float s_input[3*WARP*WARP];

    float M, S, j;

    Sum_partials_nonregular( &M, &S, &j, d_input, s_input, size);

    if(d_pp[2]>0){
      Merge(&M, &S, &j, d_pp[0], d_pp[1], d_pp[2]);
    }

    //----------------------------------------------
    //---- Writing data
    if (threadIdx.x == 0) {
      d_MSD[0] = M / j;
      d_MSD[1] = sqrt(S / j);
      d_MSD[2] = j;
      d_pp[0] = M;
      d_pp[1] = S;
      d_pp[2] = j;
    }
  }

  __global__ void MSD_GPU_Interpolate_linear(float *d_MSD_DIT, float *d_MSD_interpolated, int *d_MSD_DIT_widths, int MSD_DIT_size, int *boxcar, int max_width_performed){

    int tid  = threadIdx.x;
    if(boxcar[tid] <= max_width_performed) {
      //      int f = threadIdx.x;
      int desired_width = boxcar[tid];
      int position = (int) floorf(log2f((float) desired_width));

      float width1 = d_MSD_DIT_widths[position];
      float mean1 = d_MSD_DIT[(position)*MSD_RESULTS_SIZE];
      float StDev1 = d_MSD_DIT[(position)*MSD_RESULTS_SIZE +1];

      //      printf("\nBoxcar: %f \t desired: %f", (float)boxcar[f], desired_width);

      if(position == MSD_DIT_size-1 && width1==(int) desired_width) {
	//                (*mean) = mean1;
	//                (*StDev) = StDev1;
	d_MSD_interpolated[tid*2] = mean1;
	d_MSD_interpolated[tid*2+1] = StDev1;
      }
      else {
	float width2 = d_MSD_DIT_widths[position+1];
	float distance_in_width = width2 - width1;

	float mean2 = d_MSD_DIT[(position+1)*MSD_RESULTS_SIZE];
	float distance_in_mean = mean2 - mean1;

	float StDev2 = d_MSD_DIT[(position+1)*MSD_RESULTS_SIZE +1];
	float distance_in_StDev = StDev2 - StDev1;

	//                        printf("Position: \t %i \t f: %i\n", position, f);
	//                        printf("width:[%f;%f]; mean:[%f;%f]; sd:[%f;%f]\n",width1, width2, mean1, mean2, StDev1, StDev2);
	//                        printf("d width %f; d mean: %f; d StDef: %f\n", distance_in_width, distance_in_mean, distance_in_StDev);
	//                        printf("\tDesired_width: %f\n", desired_width);

	//                (*mean) = mean1 + (distance_in_mean/distance_in_width)*((float) desired_width - width1);
	//                (*StDev) = StDev1 + (distance_in_StDev/distance_in_width)*((float) desired_width - width1);
	d_MSD_interpolated[tid*2] = mean1 + (distance_in_mean/distance_in_width)*((float) desired_width - width1);
	d_MSD_interpolated[tid*2+1] = StDev1 + (distance_in_StDev/distance_in_width)*((float) desired_width - width1);

      }
    }
  }

  /** \brief Kernel wrapper function to MSD_GPU_final_regular kernel function. */
  void call_kernel_MSD_GPU_final_regular(const dim3 &grid_size, const dim3 &block_size, float *const d_input, float *const d_output, const int &size) {
    MSD_GPU_final_regular<<<grid_size, block_size>>>(d_input, d_output, size);
  }

  /** \brief Kernel wrapper function to MSD_GPU_final_regular kernel function. */
  void call_kernel_MSD_GPU_final_regular(const dim3 &grid_size, const dim3 &block_size, float *const d_input, float *const d_MSD, float *const d_pp, const int &size) {
    MSD_GPU_final_regular<<<grid_size,block_size>>>(d_input, d_MSD, d_pp, size);
  }

  /** \brief Kernel wrapper function to MSD_GPU_final_nonregular kernel function. */
  void call_kernel_MSD_GPU_final_nonregular(const dim3 &grid_size, const dim3 &block_size, float *const d_input, float *const d_MSD, const int &size) {
    MSD_GPU_final_nonregular<<<grid_size, block_size>>>(d_input, d_MSD, size);
  }

  /** \brief Kernel wrapper function to MSD_GPU_final_nonregular kernel function. */
  void call_kernel_MSD_GPU_final_nonregular(const dim3 &grid_size, const dim3 &block_size, float *const d_input, float *const d_MSD, float *const d_pp, const int &size) {
    MSD_GPU_final_nonregular<<<grid_size, block_size>>>(d_input, d_MSD, d_pp, size);
  }

  /** \brief Kernel wrapper function to MSD_GPU_Interpolate_linear kernel function. */
  void call_kernel_MSD_GPU_Interpolate_linear(const dim3 &grid_size, const dim3 &block_size, float *const d_MSD_DIT, float *const d_MSD_interpolated, int *const d_MSD_DIT_widths, const int &MSD_DIT_size, int *const boxcar, const int &max_width_performed) {
    MSD_GPU_Interpolate_linear<<<grid_size, block_size>>>(d_MSD_DIT, d_MSD_interpolated, d_MSD_DIT_widths, MSD_DIT_size, boxcar, max_width_performed);
  }

} //namespace astroaccelerate
