#ifndef ASTRO_ACCELERATE_AA_DEVICE_SAVE_DATA_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_SAVE_DATA_HPP

namespace astroaccelerate {
  /** \brief Save a chunk of data of the provided size from the device to the host. */
  void save_data(float *device_pointer, float *host_pointer, size_t size);

  /** \brief Save a chunk of data of the provided size from the device to the host at a given offset. */
  void save_data_offset(float *device_pointer, size_t device_offset, float *host_pointer, size_t host_offset, size_t size);

  /** \brief Asynchronous save a chunk of data of the provided size from the device to the host at a given offset. */
  void save_data_offset_stream(int dm_range, int current_time_chunk, int **t_processed, long int inc, int *inBin, int const*const ndms, float *d_DDTR_output, float ***output_buffer);
} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_SAVE_DATA_HPP

