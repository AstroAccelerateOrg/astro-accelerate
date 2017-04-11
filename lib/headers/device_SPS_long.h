#ifndef __SPS_LONG__
#define __SPS_LONG__

size_t Get_memory_requirement_of_SPS();
extern void PD_SEARCH_LONG_init();
extern int PD_SEARCH_LONG(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, int *d_output_taps, float *d_MSD, int max_boxcar_width, int nDMs, int nTimesamples, int *t_max_iterarion);
extern int PD_SEARCH_LONG_BLN(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, int max_boxcar_width, int nDMs, int nTimesamples, int *t_max_iterarion);
extern int PD_SEARCH_LONG_BLN_IF(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, int max_boxcar_width, int nDMs, int nTimesamples, int *t_max_iterarion);
#endif
