#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "aa_host_utilities.hpp"

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



} //namespace astroaccelerate