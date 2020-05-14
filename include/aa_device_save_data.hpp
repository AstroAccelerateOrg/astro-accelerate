#ifndef ASTRO_ACCELERATE_AA_DEVICE_SAVE_DATA_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_SAVE_DATA_HPP

namespace astroaccelerate {
  /** \brief Save a chunk of data of the provided size from the device to the host. */
  void save_data(float *device_pointer, float *host_pointer, size_t size);

  /** \brief Save a chunk of data of the provided size from the device to the host at a given offset. */
  void save_data_offset(float *device_pointer, size_t device_offset, float *host_pointer, size_t host_offset, size_t size);
} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_SAVE_DATA_HPP

