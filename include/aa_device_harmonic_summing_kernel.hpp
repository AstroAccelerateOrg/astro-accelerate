#ifndef ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_KERNEL_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_KERNEL_HPP

namespace astroaccelerate {
  class HRMS_ConstParams {
  public:
    static const int warp = 32;
    static const int nThreads = 256;
  };
  
  class HRMS_normal : public HRMS_ConstParams {
  public:
    static const bool remove_scalloping_loss = false;
  };
  
  class HRMS_remove_scalloping_loss : public HRMS_ConstParams {
  public:
    static const bool remove_scalloping_loss = true;
  };

  /** \brief Kernel wrapper function for PHS_GPU_kernel kernel function. */
  void call_kernel_simple_harmonic_sum_GPU_kernel(
      const dim3 &grid_size,
      const dim3 &block_size,
      float const *const d_input,
      float *const d_output_SNR,
      ushort *const d_output_harmonics,
      float *const d_MSD,
      const int &nTimesamples,
      const int &nSpectra,
      const int &nHarmonics
  );
  
  void call_kernel_greedy_harmonic_sum_GPU_kernel(
      const dim3 &grid_size,
      const dim3 &block_size,
      float const *const d_input,
      float *const d_output_SNR,
      ushort *const d_output_harmonics,
      float *const d_MSD,
      const int &nTimesamples,
      const int &nSpectra,
      const int &nHarmonics,
      bool enable_scalloping_loss_removal
  );
  
  void call_kernel_presto_plus_harmonic_sum_GPU_kernel(
      const dim3 &grid_size,
      const dim3 &block_size,
      float const *const d_input,
      float *const d_output_SNR,
      ushort *const d_output_harmonics,
      float *const d_MSD,
      const int &nTimesamples,
      const int &nSpectra,
      const int &nHarmonics,
      bool enable_scalloping_loss_removal
  );
  
  void call_kernel_presto_harmonic_sum_GPU_kernel(
      const dim3 &grid_size,
      const dim3 &block_size,
      float const *const d_input,
      float *const d_output_SNR,
      ushort *const d_output_harmonics,
      float *const d_MSD,
      const int &nTimesamples,
      const int &nSpectra,
      const int &nHarmonicsFactor,
      bool enable_scalloping_loss_removal
  );

  void call_two_dimensional_greedy_harmonic_sum_GPU_kernel(
      const dim3 &grid_size,
      const dim3 &block_size,
      float const *const d_input,
      float *const d_output_max, 
      float *const d_output_SNR,
      ushort *const d_output_harmonics,
      float *const d_MSD,
      size_t const &N_f,
      size_t const &N_fdot,
      size_t const &max_f_idx,
      size_t const &max_fdot_idx,
      size_t const &min_fdot_idx,
      size_t const nHarmonics
    );

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_HARMONIC_SUMMING_KERNEL_HPP
