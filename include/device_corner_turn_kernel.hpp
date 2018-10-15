#ifndef ASTRO_ACCELERATE_DEVICE_CORNER_TURN_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_CORNER_TURN_KERNEL_HPP

void call_kernel_simple_corner_turn_kernel(const dim3 &          block_size,
                                           const dim3 &          grid_size,
                                           unsigned short *const d_input,
                                           float *const          d_output,
                                           const int &           primary_size,
                                           const int &secondary_size);
void call_kernel_simple_corner_turn_kernel(const dim3 & block_size,
                                           const dim3 & grid_size,
                                           float *const d_input,
                                           float *const d_output,
                                           const int &  primary_size,
                                           const int &  secondary_size);
void call_kernel_corner_turn_SM_kernel(const dim3 &       grid_size,
                                       const dim3 &       block_size,
                                       float const *const d_input,
                                       float *const       d_output,
                                       const int &        primary_size,
                                       const int &        secondary_size);
void call_kernel_swap(const dim3 &          block_size,
                      const dim3 &          grid_size,
                      unsigned short *const d_input,
                      float *const          d_output,
                      const int &           nchans,
                      const int &           nsamp);
#endif
