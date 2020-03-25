#ifndef ASTRO_ACCELERATE_AA_DEVICE_ANALYSIS_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_ANALYSIS_HPP

namespace astroaccelerate {

  /**
   * \struct analysis_pulse
   * \brief Struct that contains a single analysis pulse output.
   * \details The data contains vector data referring to each peak found:
   * \details 1. Time in seconds since beginning of fil file. [seconds].
   * \details 2. Dispersion measure [parsec / cm^3].
   * \details 3. Signal/Noise ratio [dimensionless].
   * \details 4. Pulse width in number of samples [dimensionless]. 
   */
  struct analysis_pulse {
    float dispersion_measure;
    float time;
    float snr;
    float pulse_width;
  };

  /**
   * \struct analysis_output
   * \brief Struct that contains the pulses found by the analysis component.
   * \details The pulses vector contains data referring to each pulse / peak found.
   */
  struct analysis_output {
    std::vector<analysis_pulse> pulses;
    float dm_low;
    float dm_high;
  };

  /** \brief Function that performs analysis component on the GPU. */  
  bool analysis_GPU(unsigned int *h_peak_list_DM, unsigned int *h_peak_list_TS, float *h_peak_list_SNR, unsigned int *h_peak_list_BW, size_t *peak_pos, size_t max_peak_size, int i, float tstart, int t_processed, int inBin, int *maxshift, int max_ndms, int const*const ndms, float cutoff, float sigma_constant, float max_boxcar_width_in_sec, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp, int candidate_algorithm, float *d_MSD_workarea, unsigned short *d_output_taps, float *d_MSD_interpolated, unsigned long int maxTimeSamples, int enable_msd_baselinenoise, const bool dump_to_disk, const bool dump_to_user, analysis_output &output);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_DEVICE_ANALYSIS_HPP





