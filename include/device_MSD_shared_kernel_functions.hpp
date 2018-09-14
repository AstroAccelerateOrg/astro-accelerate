#ifndef ASTRO_ACCELERATE_DEVICE_MSD_SHARED_KERNEL_FUNCTIONS_HPP
#define ASTRO_ACCELERATE_DEVICE_MSD_SHARED_KERNEL_FUNCTIONS_HPP

void call_kernel_MSD_GPU_final_regular(dim3 grid_size, dim3 block_size, float *d_input, float *d_output, int size);
void call_kernel_MSD_GPU_final_regular(dim3 grid_size, dim3 block_size, float *d_input, float *d_MSD, float *d_pp, int size);
void call_kernel_MSD_GPU_final_nonregular(dim3 grid_size, dim3 block_size, float *d_input, float *d_MSD, int size);
void call_kernel_MSD_GPU_final_nonregular(dim3 grid_size, dim3 block_size, float *d_input, float *d_MSD, float *d_pp, int size);
void call_kernel_MSD_GPU_Interpolate_linear(dim3 grid_size, dim3 block_size, float *d_MSD_DIT, float *d_MSD_interpolated, int *d_MSD_DIT_widths, int MSD_DIT_size, int *boxcar, int max_width_performed);
#endif
