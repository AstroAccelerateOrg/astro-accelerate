#ifndef ASTRO_ACCELERATE_DEVICE_RFI_KERNEL_HPP
#define ASTRO_ACCELERATE_DEVICE_RFI_KERNEL_HPP

void call_kernel_rfi_gpu_kernel(dim3 grid_size, dim3 block_size,
				unsigned short *d_input,
				int nchans, int nsamp);

#endif
