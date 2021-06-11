#ifndef ASTRO_ACCELERATE_AA_HOST_UTILITIES_HPP
#define ASTRO_ACCELERATE_AA_HOST_UTILITIES_HPP

#include <iostream>
#include <cufft.h>
#include <stdio.h>

namespace astroaccelerate {




void MSD_Kahan(float *h_input, size_t nDMs, size_t nTimesamples, size_t offset, double *mean, double *sd);

float Calculate_median(float *data, size_t primary_size);

void Export_data_to_file(float2 *data, size_t primary_size, size_t secondary_size, const char *filename);

void Export_data_to_file(float *data, size_t primary_size, size_t secondary_size, const char *filename);

void Export_data_to_file(float2 *data1, float2 *data2, float2 *data3, size_t primary_size, const char *filename);

void dered_with_MSD(float2 *data, int nSamples, int *segment_sizes, int nSegments, float *MSD);

void dered_with_MED(float2 *data, int nSamples, int *segment_sizes, int nSegments, float *MED);

void CPU_spectral_whitening(float* d_FFT_complex_output, float *d_dedispersed_data, float t_dm_step, float t_dm_low, size_t t_nTimesamples, size_t t_nTSamplesFFT, size_t t_nDMs_per_batch, size_t t_DM_shift);


} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_HOST_UTILITIES_HPP