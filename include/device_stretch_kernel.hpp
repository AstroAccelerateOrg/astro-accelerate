#ifndef ASTRO_ACCELERATE_DEVICE_STRETCH_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_STRETCH_KERNEL_HPP

void call_kernel_stretch_kernel(const dim3 &        block_size,
                                const dim3 &        grid_size,
                                const int &         smem_bytes,
                                const cudaStream_t &stream,
                                const int &         acc,
                                const int &         samps,
                                const float &       tsamp,
                                float *const        d_input,
                                float *const        d_output,
                                const float &       t_zero,
                                const float &       multiplier,
                                const float &       tsamp_inverse);

#endif
