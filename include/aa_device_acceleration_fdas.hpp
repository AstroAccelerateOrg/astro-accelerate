#ifndef ASTRO_ACCELERATE_AA_DEVICE_ACCELERATION_FDAS_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_ACCELERATION_FDAS_HPP

namespace astroaccelerate {

  /** \brief Function that performs fourier domain accelerated search (fdas). */
  void acceleration_fdas(int range,
			 int nsamp,
			 int max_ndms,
			 int processed,
			 float cutoff,
			 float ***output_buffer,
			 int const*const ndms,
			 int *inBin,
			 float *dm_low,
			 float *dm_high,
			 float *dm_step,
			 float tsamp,
			 const bool enable_custom_fft,
			 const bool enable_inbin,
			 const bool enable_norm,
			 float sigma_constant,
			 const bool enable_output_ffdot_plan,
			 const bool enable_output_fdas_list);

} // namespace astroaccelerate
#endif // ASTRO_ACCELERATE_AA_DEVICE_ACCELERATION_FDAS_HPP
