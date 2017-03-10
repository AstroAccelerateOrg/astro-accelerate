/* FDAS device functions header */
#ifndef ASTROACCELERATE_FDAS_DEVICE_H_
#define ASTROACCELERATE_FDAS_DEVICE_H_

#include <cufft.h>

/* CUDA kernels */

__global__ void cuda_overlap_copy(float2* d_ext_data, float2* d_cpx_signal, int sigblock,  int sig_rfftlen, int sig_tot_convlen, int kern_offset, int total_blocks);

__global__ void cuda_overlap_copy_smallblk(float2* d_ext_data, float2* d_cpx_signal, int sigblock,  int sig_rfftlen, int sig_tot_convlen, int kern_offset, int total_blocks);

__global__ void cuda_convolve_reg_1d_halftemps(float2* d_kernel, float2* d_signal, float2* d_ffdot_plane, int sig_tot_convlen, float scale);

__global__ void cuda_ffdotpow_concat_2d(float2* d_ffdot_plane_cpx, float* d_ffdot_plane, int sigblock, int kern_offset, int total_blocks,  int sig_tot_convlen, int sig_totlen);

__global__ void cuda_ffdotpow_concat_2d_inbin(float2* d_ffdot_plane_cpx, float* d_ffdot_plane, int sigblock, int kern_offset, int total_blocks, int  sig_tot_convlen, int sig_totlen);

__global__ void cuda_ffdotpow_concat_2d_ndm_inbin(float2* d_ffdot_plane_cpx, float* d_ffdot_plane, int kernlen, int siglen, int nkern, int kern_offset, int total_blocks, int sig_tot_convlen, int sig_totlen, int ndm);

__global__ void cuda_convolve_customfft_wes_no_reorder02(float2* d_kernel, float2* d_signal, float *d_ffdot_pw, int sigblock, int sig_tot_convlen, int sig_totlen, int offset, float scale);

__global__ void cuda_convolve_customfft_wes_no_reorder02_inbin(float2* d_kernel, float2* d_signal, float *d_ffdot_pw, int sigblock, int sig_tot_convlen, int sig_totlen, int offset, float scale, float2 *ip_edge_points);

__global__ void customfft_fwd_temps_no_reorder(float2* d_signal);

#endif /* FDAS_DEVICE_H */
