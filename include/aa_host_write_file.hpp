#ifndef ASTRO_ACCELERATE_AA_HOST_WRITE_FILE_HPP
#define ASTRO_ACCELERATE_AA_HOST_WRITE_FILE_HPP

namespace astroaccelerate {

  /** \brief This function writes the transformed space (dm,t) out to file in binary format. */
  void write_output(int i, int t_processed, int ndms, float *output_buffer, float *dm_low, float *dm_high);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_HOST_WRITE_FILE_HPP
