#ifndef ASTRO_ACCELERATE_AA_DEVICE_ANALYSIS_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_ANALYSIS_HPP

namespace astroaccelerate {

  struct analysis_output {
    /**
     * The std::vector<float> contains 4 elements per processed dm range:
     * 1. Time in seconds since beginning of fil file. [seconds].
     * 2. Dispersion measure [parsec / cm^3].
     * 3. Signal/Noise ratio [dimensionless].
     * 4. Pulse width in number of samples [dimensionless].
     */
    std::vector<float> data;
    float              dm_low;
    float              dm_high;
  };
  
  void analysis_GPU(unsigned int *h_peak_list_DM, unsigned int *h_peak_list_TS, float *h_peak_list_SNR, unsigned int *h_peak_list_BW, size_t *peak_pos, size_t max_peak_size, int i, float tstart, int t_processed, int inBin, int *maxshift, int max_ndms, int const*const ndms, float cutoff, float OR_sigma_multiplier, float max_boxcar_width_in_sec, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp, int candidate_algorithm, float *d_MSD_workarea, unsigned short *d_output_taps, float *d_MSD_interpolated, unsigned long int maxTimeSamples, int enable_sps_baselinenoise, const bool dump_to_disk, const bool dump_to_user, analysis_output &output);

} //namespace astroaccelerate
  
#endif /* ASTRO_ACCELERATE_AA_DEVICE_ANALYSIS_HPP */





