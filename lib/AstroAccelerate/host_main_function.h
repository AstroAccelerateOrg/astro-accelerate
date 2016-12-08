#ifndef SKA_ASTROACCELERATE_MAIN_FUNCTION_H_
#define SKA_ASTROACCELERATE_MAIN_FUNCTION_H_

void main_function
(
	int argc,
	char* argv[],
	// Internal code variables
	// Counters and flags
	int range,
	int enable_debug,
	int enable_analysis,
	int enable_acceleration,
	int enable_periodicity,
	int output_dmt,
	int enable_zero_dm,
	int *inBin,
	int *outBin,
	int *ndms,
	int maxshift,
	int max_ndms,
	int max_samps,
	int num_tchunks,
	int total_ndms,
	int multi_file,
	float max_dm,
	// Memory sizes and pointers
    size_t inputsize,
    size_t outputsize,
	size_t gpu_inputsize,
	size_t gpu_outputsize,
    unsigned short  *input_buffer,
	float ***output_buffer,
	unsigned short  *d_input,
	float *d_output,
	float *dmshifts,
	float *user_dm_low,
	float *user_dm_high,
	float *user_dm_step,
	float *dm_low,
	float *dm_high,
	float *dm_step,
	// Telescope parameters
	int nchans,
	int nsamp,
	int nbits,
	int nsamples,
	int nifs,
	int **t_processed,
	int nboots,
	int ntrial_bins,
	int navdms,
	int nsearch,
	float aggression,
	float narrow,
	float wide,
	int	maxshift_original,
	double	tsamp_original,
	float tstart,
	float tsamp,
	float fch1,
	float foff,
	// Analysis variables
	float power,
	float sigma_cutoff,
	double start_time
);

#endif
