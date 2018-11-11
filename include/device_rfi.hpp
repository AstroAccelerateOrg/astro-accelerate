#ifndef ASTRO_ACCELERATE_DEVICE_RFI_HPP
#define ASTRO_ACCELERATE_DEVICE_RFI_HPP

namespace astroaccelerate {

extern void rfi_gpu(unsigned short *d_input, int nchans, int nsamp);

} //namespace astroaccelerate
  
#endif

