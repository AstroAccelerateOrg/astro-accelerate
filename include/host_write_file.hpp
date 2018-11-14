#ifndef ASTRO_ACCELERATE_HOST_WRITE_FILE_HPP
#define ASTRO_ACCELERATE_HOST_WRITE_FILE_HPP

namespace astroaccelerate {

void write_output(int i, int t_processed, int ndms, float *output_buffer, float *dm_low, float *dm_high);

} //namespace astroaccelerate
  
#endif
