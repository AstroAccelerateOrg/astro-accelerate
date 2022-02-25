#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cuda.h>
#include <cuda_runtime.h>
#include "aa_host_utilities.hpp"
#include "aa_device_spectrum_whitening.hpp"
#include "presto_funcs.hpp"
#include "presto.hpp"

namespace astroaccelerate {

//---------------------------------------------------------------------------------
//-------> Kahan MSD
void d_kahan_summation(float *signal, size_t nDMs, size_t nTimesamples, size_t offset, float *result, float *error){
	double sum;
	double sum_error;
	double a,b;
	
	sum=0;
	sum_error=0;
	for(size_t d=0;d<nDMs; d++){
		for(size_t s=0; s<(nTimesamples-offset); s++){
			a=signal[(size_t) (d*nTimesamples + s)]-sum_error;
			b=sum+a;
			sum_error=(b-sum);
			sum_error=sum_error-a;
			sum=b;
		}
	}
	*result=sum;
	*error=sum_error;
}

void d_kahan_sd(float *signal, size_t nDMs, size_t nTimesamples, size_t offset, double mean, float *result, float *error){
	double sum;
	double sum_error;
	double a,b,dtemp;
	
	sum=0;
	sum_error=0;
	for(size_t d=0;d<nDMs; d++){
		for(size_t s=0; s<(nTimesamples-offset); s++){
			dtemp=(signal[(size_t) (d*nTimesamples + s)]-sum_error - mean);
			a=dtemp*dtemp;
			b=sum+a;
			sum_error=(b-sum);
			sum_error=sum_error-a;
			sum=b;
		}
	}
	*result=sum;
	*error=sum_error;
}

void MSD_Kahan(float *h_input, size_t nDMs, size_t nTimesamples, size_t offset, double *mean, double *sd){
	float error, signal_mean, signal_sd;
	size_t nElements=nDMs*(nTimesamples-offset);
	
	d_kahan_summation(h_input, nDMs, nTimesamples, offset, &signal_mean, &error);
	signal_mean=signal_mean/nElements;
	
	d_kahan_sd(h_input, nDMs, nTimesamples, offset, signal_mean, &signal_sd, &error);
	signal_sd=sqrt(signal_sd/nElements);

	*mean=signal_mean;
	*sd=signal_sd;
}
//-------> Kahan MSD
//---------------------------------------------------------------------------------

float Calculate_median(float *data, size_t primary_size){
	float median;
	std::sort(data, data + primary_size);
	int middle = (primary_size>>1);
	if(primary_size%2==0){
		median = (data[middle] + data[middle + 1])/2;
		//median = data[middle + 1];
	}
	else {
		median = data[middle];
	}
	return(median);
}


void Export_data_to_file(float2 *data, size_t primary_size, size_t secondary_size, const char *filename){
    std::ofstream FILEOUT;
    FILEOUT.open (filename, std::ofstream::out);
	for(size_t sec=0; sec<secondary_size; sec++){
		for(size_t prim=0; prim<primary_size; prim++){
			size_t pos = sec*primary_size + prim;
			double power = sqrt(data[pos].x*data[pos].x + data[pos].y*data[pos].y);
			FILEOUT << data[pos].x << " " << data[pos].y << " " << power << std::endl;
		}
		FILEOUT << std::endl;
		FILEOUT << std::endl;
	}
    FILEOUT.close();
}

void Export_data_to_file(float *data, size_t primary_size, size_t secondary_size, const char *filename){
    std::ofstream FILEOUT;
    FILEOUT.open (filename, std::ofstream::out);
	for(size_t sec=0; sec<secondary_size; sec++){
		for(size_t prim=0; prim<primary_size; prim++){
			size_t pos = sec*primary_size + prim;
			FILEOUT << data[pos] << std::endl;
		}
		FILEOUT << std::endl;
		FILEOUT << std::endl;
	}
    FILEOUT.close();
}

void Export_data_to_file(float2 *data1, float2 *data2, float2 *data3, size_t primary_size, const char *filename){
    std::ofstream FILEOUT;
    FILEOUT.open (filename, std::ofstream::out);
	for(size_t prim=0; prim<primary_size; prim++){
		size_t pos = prim;
		double power1 = sqrt(data1[pos].x*data1[pos].x + data1[pos].y*data1[pos].y);
		double power2 = sqrt(data2[pos].x*data2[pos].x + data2[pos].y*data2[pos].y);
		double power3 = sqrt(data3[pos].x*data3[pos].x + data3[pos].y*data3[pos].y);
		FILEOUT << power1 << " " << power2 << " " << power3 << std::endl;
	}
	FILEOUT << std::endl;
	FILEOUT << std::endl;
    FILEOUT.close();
}

void dered_with_MSD(float2 *data, int nSamples, int *segment_sizes, int nSegments, float *MSD){
	for(int s=0; s<nSegments; s++){
		if(s==0){
			// First segment cannot calculate slope thus it is normalized by mean only
			int f0_pos  = segment_sizes[s + 0];
			int fp1_pos = segment_sizes[s + 1];
			int current_range  = fp1_pos - f0_pos;
			float current_mean = MSD[2*s];
			float norm = 1.0/sqrt(current_mean);
			
			data[0].x = 1.0;
			data[0].y = 0.0;
			for(int i=1; i<(current_range>>1); i++){
				data[i].x *= norm;
				data[i].y *= norm;
			}
		}
		else {
			// This is presto deredning scheme
			int fm1_pos = segment_sizes[s - 1];
			int f0_pos  = segment_sizes[s + 0];
			int fp1_pos = segment_sizes[s + 1];
			int previous_range = f0_pos - fm1_pos;
			int current_range  = fp1_pos - f0_pos;
			float previous_mean = MSD[2*(s - 1)];
			float current_mean  = MSD[2*(s)];
			int range = ((previous_range + current_range)>>1);
			float slope = (current_mean - previous_mean) / ((float) range);
			
			for(int i=0; i<=range; i++){
				int local_pos = fm1_pos + (previous_range>>1) + i;
				float norm = 1.0/sqrt(previous_mean + slope*(range-i));
				if(local_pos < nSamples){
					data[local_pos].x *= norm;
					data[local_pos].y *= norm;
				}
			}
			
			// deal with the end
			if( s==(nSegments-2) ){
				int range_to_end = nSamples - (fm1_pos + (previous_range>>1) + range);
				for(int i=0; i<range_to_end; i++){
					int local_pos = fm1_pos + (previous_range>>1) + range + i;
					float norm = 1.0/sqrt(previous_mean + slope*i);
					if(local_pos < nSamples){
						data[local_pos].x *= norm;
						data[local_pos].y *= norm;
					}
				}
			}
		}
	}
}

void dered_with_MED(float2 *data, int nSamples, int *segment_sizes, int nSegments, float *MED){
	for(int s=0; s<nSegments; s++){
		if(s==0){
			// First segment cannot calculate slope thus it is normalized by mean only
			int f0_pos  = segment_sizes[s + 0];
			int fp1_pos = segment_sizes[s + 1];
			int current_range  = fp1_pos - f0_pos;
			float current_mean = MED[s] / log(2.0);
			float norm = 1.0/sqrt(current_mean);
			
			data[0].x = 1.0;
			data[0].y = 0.0;
			for(int i=1; i<(current_range>>1); i++){
				data[i].x *= norm;
				data[i].y *= norm;
			}
		}
		else {
			// This is presto deredning scheme
			int fm1_pos = segment_sizes[s - 1];
			int f0_pos  = segment_sizes[s + 0];
			int fp1_pos = segment_sizes[s + 1];
			int previous_range = f0_pos - fm1_pos;
			int current_range  = fp1_pos - f0_pos;
			float previous_mean = MED[s - 1] / log(2.0);
			float current_mean  = MED[s] / log(2.0);
			int range   = ((previous_range + current_range)>>1);
			float slope = (current_mean - previous_mean) / ((float) range);
			//printf("previous_mean=%e; current_mean=%e; slope=%e;\n", previous_mean, current_mean, slope);
			for(int i=0; i<=range; i++){
				int local_pos = fm1_pos + (previous_range>>1) + i;
				float norm = 1.0/sqrt(previous_mean + slope*(range-i));
				if(local_pos < nSamples){
					data[local_pos].x *= norm;
					data[local_pos].y *= norm;
				}
			}
			
			// deal with the end
			if( s==(nSegments-2) ){
				int range_to_end = nSamples - (fm1_pos + (previous_range>>1) + range);
				for(int i=0; i<range_to_end; i++){
					int local_pos = fm1_pos + (previous_range>>1) + range + i;
					float norm = 1.0/sqrt(previous_mean + slope*i);
					if(local_pos < nSamples){
						data[local_pos].x *= norm;
						data[local_pos].y *= norm;
					}
				}
			}
		}
	}
}

void CPU_spectral_whitening(float* d_FFT_complex_output, float *d_dedispersed_data, float t_dm_step, float t_dm_low, size_t t_nTimesamples, size_t t_nTSamplesFFT, size_t t_nDMs_per_batch, size_t t_DM_shift){
        
    //---------> CPU spectral whitening
    printf("full nTimesamples=%zu; half nTimesamples=%zu;\n", t_nTimesamples, t_nTSamplesFFT);
    // Copy stuff to the host
    cudaError_t err;
    char filename[300];
    size_t fft_input_size_bytes = t_nTSamplesFFT*t_nDMs_per_batch*sizeof(float2);
    size_t fft_power_size_bytes = t_nTSamplesFFT*t_nDMs_per_batch*sizeof(float);
    size_t fft_ddtr_size_bytes  = t_nTimesamples*t_nDMs_per_batch*sizeof(float);
    float2 *h_fft_input;
    float  *h_fft_power;
    float  *h_ddtr_data;
    printf("Data copied and power calculation...\n");
    h_fft_input = (float2*) malloc(fft_input_size_bytes);
    h_fft_power = (float*) malloc(fft_power_size_bytes);
    h_ddtr_data = (float*) malloc(fft_ddtr_size_bytes);
    
    err = cudaMemcpy(h_fft_input, d_FFT_complex_output, fft_input_size_bytes, cudaMemcpyDeviceToHost);
    if(err != cudaSuccess) printf("CUDA error\n");
    err = cudaMemcpy(h_ddtr_data, d_dedispersed_data, fft_ddtr_size_bytes, cudaMemcpyDeviceToHost);
    if(err != cudaSuccess) printf("CUDA error\n");
    
    for(size_t d=0; d<t_nDMs_per_batch; d++){
        for(size_t s=0; s<t_nTSamplesFFT; s++){
            size_t pos = d*t_nTSamplesFFT + s;
            h_fft_power[pos] = h_fft_input[pos].x*h_fft_input[pos].x + h_fft_input[pos].y*h_fft_input[pos].y;
        }
    }
    
    
    // Export data to file
    printf("Exporting fft data to file...\n");
    for(size_t d=0; d<t_nDMs_per_batch; d++){
        sprintf(filename, "PSR_fft_data_%f.dat", t_dm_low + t_dm_step*(t_DM_shift + d));
        size_t pos = d*t_nTSamplesFFT;
        Export_data_to_file(&h_fft_input[pos], t_nTSamplesFFT, 1, filename);
    }
    printf("Exporting ddtr data to file...\n");
    for(size_t d=0; d<t_nDMs_per_batch; d++){
        sprintf(filename, "PSR_ddtr_data_%f.dat", t_dm_low + t_dm_step*(t_DM_shift + d));
        size_t pos = d*t_nTimesamples;
        Export_data_to_file(&h_ddtr_data[pos], t_nTimesamples, 1, filename);
    }
    
    
    // Create segments for de-redning
    printf("Calculating segment sizes...\n");
    int max_segment_length = 256;
    int min_segment_length = 6;
    std::vector<int> segment_sizes;
    create_dered_segment_sizes_prefix_sum(&segment_sizes, min_segment_length, max_segment_length, t_nTSamplesFFT);
    int nSegments = segment_sizes.size();
    
    // Calculate mean for segments and export_size
    printf("Calculating MSD...\n");
    size_t MSD_segmented_size_bytes = 2*nSegments*t_nDMs_per_batch*sizeof(float);
    float *h_segmented_MSD;
    h_segmented_MSD = (float*) malloc(MSD_segmented_size_bytes);
    for(size_t d=0; d<t_nDMs_per_batch; d++){
        for(int s=0; s<(nSegments - 1); s++){
            size_t MSD_pos = d*nSegments + s;
            double mean, stdev;
            int range = segment_sizes[s + 1] - segment_sizes[s];
            size_t pos = d*t_nTSamplesFFT + segment_sizes[s];
            MSD_Kahan(&h_fft_power[pos], 1, range, 0, &mean, &stdev);
            h_segmented_MSD[2*MSD_pos] = (float) mean;
            h_segmented_MSD[2*MSD_pos + 1] = (float) stdev;
        }
        
        sprintf(filename, "PSR_fft_data_means_%f.dat", t_dm_low + t_dm_step*(t_DM_shift + d));
        Export_data_to_file((float2*) &h_segmented_MSD[d*nSegments], nSegments, 1, filename);
    }
    
    // Calculate medians for segments
    printf("Calculating median...\n");
    size_t MED_segmented_size_bytes = nSegments*t_nDMs_per_batch*sizeof(float);
    float *h_segmented_MED;
    h_segmented_MED = (float*) malloc(MED_segmented_size_bytes);
    for(size_t d=0; d<t_nDMs_per_batch; d++){
        for(int s=0; s<(nSegments - 1); s++){
            size_t MED_pos = d*nSegments + s;
            int range = segment_sizes[s + 1] - segment_sizes[s];
            size_t pos = d*t_nTSamplesFFT + segment_sizes[s];
            h_segmented_MED[MED_pos] = Calculate_median(&h_fft_power[pos], range);
        }
        
        sprintf(filename, "PSR_fft_data_median_%f.dat", t_dm_low + t_dm_step*(t_DM_shift + d));
        Export_data_to_file(&h_segmented_MED[d*nSegments], nSegments, 1, filename);
    }
    
    // Calculate medians for segments based on presto
    printf("Calculating median using presto...\n");
    float *h_segmented_MED_p;
    h_segmented_MED_p = (float*) malloc(MED_segmented_size_bytes);
    for(size_t d=0; d<t_nDMs_per_batch; d++){
        for(int s=0; s<(nSegments - 1); s++){
            size_t MED_pos = d*nSegments + s;
            int range = segment_sizes[s + 1] - segment_sizes[s];
            size_t pos = d*t_nTSamplesFFT + segment_sizes[s];
            h_segmented_MED_p[MED_pos] = median(&h_fft_power[pos], range);
        }
        
        sprintf(filename, "PSR_fft_data_median_p_%f.dat", t_dm_low + t_dm_step*(t_DM_shift + d));
        Export_data_to_file(&h_segmented_MED_p[d*nSegments], nSegments, 1, filename);
    }
    
    // Deredning by presto
    printf("De-redning...\n");
    for(size_t d=0; d<t_nDMs_per_batch; d++){
        float2 *presto_dered, *MSD_dered, *MED_dered;
        presto_dered = new float2[t_nTSamplesFFT];
        MSD_dered    = new float2[t_nTSamplesFFT];
        MED_dered    = new float2[t_nTSamplesFFT];
        for(size_t s=0; s<t_nTSamplesFFT; s++){
            size_t pos = d*t_nTSamplesFFT + s;
            presto_dered[s] = h_fft_input[pos];
            MSD_dered[s]    = h_fft_input[pos];
            MED_dered[s]    = h_fft_input[pos];
        }
        size_t MSD_pos = d*nSegments;
        presto_dered_sig(presto_dered, t_nTSamplesFFT);
        dered_with_MSD(MSD_dered, t_nTSamplesFFT, segment_sizes.data(), nSegments, (float*) &h_segmented_MSD[MSD_pos]);
        dered_with_MED(MED_dered, t_nTSamplesFFT, segment_sizes.data(), nSegments, &h_segmented_MED[MSD_pos]);
        sprintf(filename, "PSR_fft_dered_%f.dat", t_dm_low + t_dm_step*(t_DM_shift + d));
        Export_data_to_file(presto_dered, MSD_dered, MED_dered, t_nTSamplesFFT, filename);
    }
    
    free(h_fft_input);
    free(h_fft_power);
    free(h_segmented_MSD);
    free(h_segmented_MED);
}


} //namespace astroaccelerate