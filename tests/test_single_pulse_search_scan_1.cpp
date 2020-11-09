#include <iostream>
#include <sstream>

#include <math.h>
#include "aa_params.hpp"
#include "aa_device_info.hpp"
#include "aa_device_SPS_long_kernel.hpp"

using namespace astroaccelerate;

double  max_error = 10.e-4;


float get_error(float A, float B){
	float error, div_error=10000, per_error=10000, order=0;
	int power;
	if(A<0) A = -A;
	if(B<0) B = -B;
	
	if (A>B) {
		div_error = A-B;
		if(B>10){
			power = (int) log10(B);
			order = pow(10,power);
			div_error = div_error/order;
		}
	}
	else {
		div_error = B-A;
		if(A>10){
			power = (int) log10(A);
			order = pow(10,power);
			div_error = div_error/order;
		}
	}
	
	if(div_error<per_error) error = div_error;
	else error = per_error;
	return(error);
}

int CompareData(float *CPU_result, float *GPU_result, float CPU_scale, float GPU_scale, int CPU_offset, int GPU_offset, int CPU_dim_x, int GPU_dim_x, int dim_y, int nSamples, double *total_error, double *mean_error){
	double total_error_l = 0, mean_error_l = 0;
	size_t nErrors = 0;
	int cislo = 0;
	float error;
	printf("DEBUG: CompareData variables: CPU_dim_x=%d; GPU_dim_x=%d; dim_y=%d; nSamples=%d;\n", CPU_dim_x, GPU_dim_x, dim_y, nSamples);
	
	for(int y=0; y<dim_y; y++){
		for(int x=0; x<nSamples; x++){
			int CPU_pos = y*CPU_dim_x + x + CPU_offset;
			int GPU_pos = y*GPU_dim_x + x + GPU_offset;
			float CPU, GPU;
			CPU = CPU_result[CPU_pos]/CPU_scale;
			GPU = GPU_result[GPU_pos]/GPU_scale;
			
			error = get_error(CPU, GPU);
			total_error_l = total_error_l + error;
			if( error > max_error ){
				nErrors++;
				if(cislo<40){
					printf("Error [%f] CPU=%f; GPU=%f; x=%d; y=%d;\n", error, CPU, GPU, x, y);
					cislo++;
				}
			}
		}
	}
	mean_error_l = total_error_l/(((double) nSamples)*((double) dim_y));
	(*total_error) = total_error_l;
	(*mean_error) = mean_error_l;
	return(nErrors);
}


int CompareDataTaps(ushort *CPU_result, ushort *GPU_result, float CPU_scale, float GPU_scale, int CPU_offset, int GPU_offset, int CPU_dim_x, int GPU_dim_x, int dim_y, int nSamples, double *total_error, double *mean_error){
	double total_error_l = 0, mean_error_l = 0;
	size_t nErrors = 0;
	int cislo = 50;
	float error;
	printf("DEBUG: CompareData variables: CPU_dim_x=%d; GPU_dim_x=%d; dim_y=%d; nSamples=%d;\n", CPU_dim_x, GPU_dim_x, dim_y, nSamples);
	
	for(int y=0; y<dim_y; y++){
		for(int x=0; x<nSamples; x++){
			int CPU_pos = y*CPU_dim_x + x + CPU_offset;
			int GPU_pos = y*GPU_dim_x + x + GPU_offset;
			float CPU, GPU;
			CPU = (float) CPU_result[CPU_pos]/CPU_scale;
			GPU = (float) GPU_result[GPU_pos]/GPU_scale;
			
			error = abs(CPU - GPU);
			total_error_l = total_error_l + error;
			if( error > 1.2 ){
				nErrors++;
				if(cislo<40){
					printf("Error [%f] CPU=%f; GPU=%f; x=%d; y=%d;\n", error, CPU, GPU, x, y);
					cislo++;
				}
			}
		}
	}
	mean_error_l = total_error_l/(((double) nSamples)*((double) dim_y));
	(*total_error) = total_error_l;
	(*mean_error) = mean_error_l;
	return(nErrors);
}


int CompareData(float *CPU_result, float *GPU_result, size_t dim_x, size_t dim_y);

void Decimate_in_time(float *h_input, float *h_CPU_decimate, int DIT_value, int DIT_factor, size_t nDMs, size_t nTimesamples){
	float ftemp;
	size_t decimated_timesamples;
	
	decimated_timesamples=nTimesamples/(DIT_value*DIT_factor);
	for(size_t d=0; d<nDMs; d++){
		for(size_t s=0; s<decimated_timesamples; s++){
			ftemp=0;
			for(int t=0; t<DIT_factor; t++){
				ftemp = ftemp + h_input[d*decimated_timesamples*DIT_factor + s*DIT_factor + t];
			}
			h_CPU_decimate[d*decimated_timesamples + s]=ftemp;
		}
	}
}

void Single_pulse_search_offcentral_1st(float *h_input, float *h_MSD_interpolated, float *h_CPU_boxcar_values, float *h_output_SNR, ushort *h_output_taps, int maxTaps, size_t nDMs, size_t nTimesamples){
	int nTaps;
	size_t d, s, output_pos;
	float ftemp, SNR, FIR_value, mean, stdev;
	
	// Calculate SNR for unprocessed signal (equivalent to producing boxcar of width 1)
	for(d=0; d<nDMs; d++){
		for(s=0; s<(nTimesamples-maxTaps+1); s++){
			h_output_SNR[d*nTimesamples + s]  = (h_input[d*nTimesamples + s] - h_MSD_interpolated[0])/h_MSD_interpolated[1];
			h_output_taps[d*nTimesamples + s] = 1;
		}
	}
	
	// Calculating boxcar filters with widths up to and including max_FIR.
	size_t decimated_timesamples = nTimesamples;
	size_t dtm = decimated_timesamples>>1;
	for(d=0; d<nDMs; d++){
		for(s=0; s<(decimated_timesamples-maxTaps+1); s++){
			// loading initial data
			nTaps = 1;
			FIR_value = h_input[d*decimated_timesamples + s];
			ftemp=sqrt((float) nTaps);
			SNR=FIR_value/ftemp;
			output_pos = d*nTimesamples + s;
			if(SNR>h_output_SNR[output_pos]) {
				h_output_SNR[output_pos]=SNR;
				h_output_taps[output_pos]=nTaps;
			}
			
			for(int t=1; t<maxTaps; t++){
				nTaps = t+1;
				mean  = h_MSD_interpolated[2*t];
				stdev = h_MSD_interpolated[2*t+1];
				FIR_value = FIR_value + h_input[d*decimated_timesamples + s + t];
				
				SNR = (FIR_value-mean)/stdev;
				output_pos = d*nTimesamples + s;
				if(SNR>h_output_SNR[output_pos]) {
					h_output_SNR[output_pos]  = SNR;
					h_output_taps[output_pos] = nTaps;
				}
			}
			
			if( (s%2)==0 ) {
				h_CPU_boxcar_values[d*dtm + (int) (s/2)] = FIR_value;
			}
		}
	}
	//at this point we should have maximum SNR for each point in the DM/Timesample plane
	
	for(d=0; d<nDMs; d++){
		for(s=(decimated_timesamples-maxTaps+1); s<(decimated_timesamples); s++){
			output_pos = d*decimated_timesamples + s;
			if( (s%2)==0 && output_pos<(nDMs*decimated_timesamples) ) {
				FIR_value = h_input[d*decimated_timesamples + s];			
				for(int t=1; t<maxTaps; t++) FIR_value = FIR_value + h_input[d*decimated_timesamples + s + t];
				h_CPU_boxcar_values[d*dtm + (int) (s/2)] = FIR_value;
			}
		}
	}
	
}




int check_boxcar_results(float *h_input, float *h_MSD_interpolated, float *h_GPU_boxcar_values, float *h_GPU_decimated, float *h_GPU_output_SNR, ushort *h_GPU_output_taps, int nBoxcars, size_t nTimesamples, size_t nDMs){
	float *h_CPU_boxcar_values, *h_CPU_decimated, *h_CPU_output_SNR;
	ushort *h_CPU_output_taps;
	
	size_t dec_nTimesamples = (nTimesamples>>1);
	double total_error, mean_error;
	
	
	int decimate_pass;
	h_CPU_decimated = new float[dec_nTimesamples*nDMs];
	Decimate_in_time(h_input, h_CPU_decimated, 1, 2, nDMs, nTimesamples);
	decimate_pass = CompareData(h_CPU_decimated, h_GPU_decimated, 1.0, 1.0, 0, 0, dec_nTimesamples, dec_nTimesamples, nDMs, dec_nTimesamples-(nBoxcars>>1), &total_error, &mean_error);
	printf("DEBUG: Decimation comparison: nErrors=%d; total_error = %e; mean_error = %e;\n", decimate_pass, total_error, mean_error);
	delete[] h_CPU_decimated;
	
	int SNR_pass;
	int boxcar_value_pass;
	int taps_pass;
	h_CPU_boxcar_values = new float[dec_nTimesamples*nDMs];
	h_CPU_output_SNR    = new float[nTimesamples*nDMs*2];
	h_CPU_output_taps   = new ushort[nTimesamples*nDMs*2];
	Single_pulse_search_offcentral_1st(h_input, h_MSD_interpolated, h_CPU_boxcar_values, h_CPU_output_SNR, h_CPU_output_taps, nBoxcars, nDMs, nTimesamples);
	
	boxcar_value_pass = CompareData(h_CPU_boxcar_values, h_GPU_boxcar_values, 1.0, 1.0, 0, 0, dec_nTimesamples, dec_nTimesamples, nDMs, dec_nTimesamples-(nBoxcars>>1), &total_error, &mean_error);
	printf("DEBUG: Boxcar value comparison: nErrors=%d; total_error = %e; mean_error = %e;\n", boxcar_value_pass, total_error, mean_error);
	
	SNR_pass = CompareData(h_CPU_output_SNR, h_GPU_output_SNR, 1.0, 1.0, 0, 0, nTimesamples, nTimesamples, nDMs, nTimesamples-nBoxcars, &total_error, &mean_error);
	printf("DEBUG: Output SNR comparison: nErrors=%d; total_error = %e; mean_error = %e;\n", SNR_pass, total_error, mean_error);	

	taps_pass = CompareDataTaps(h_CPU_output_taps, h_GPU_output_taps, 1.0, 1.0, 0, 0, nTimesamples, nTimesamples, nDMs, nTimesamples-nBoxcars, &total_error, &mean_error);
	printf("DEBUG: Output taps comparison: nErrors=%d; total_error = %e; mean_error = %e;\n", taps_pass, total_error, mean_error);
	printf("DEBUG: The number of errors for the boxcar filter length might be higher then 0.\n");
	if(mean_error<1.0e-4) taps_pass = 0;
	
	
	delete[] h_CPU_boxcar_values;
	delete[] h_CPU_output_SNR;
	delete[] h_CPU_output_taps;
	
	
	if(decimate_pass==0 && boxcar_value_pass==0 && SNR_pass==0) return(1);
	else return(0);
}

int main(void) {
	cudaSetDevice(CARD);
	
	// GPU memory
	float  *d_input;
	float  *d_boxcar_values;
	float  *d_decimated;
	float  *d_output_SNR;
	ushort *d_output_taps;
	float  *d_MSD_interpolated;
	
	// Memory
	size_t free_mem,total_mem, max_mem;
	max_mem = 8l*1024l*1024l*1024l;
	cudaMemGetInfo(&free_mem,&total_mem);
	if(free_mem>(max_mem)) free_mem = max_mem; 
	int nTimesamples = 800000;
	size_t memory_required_per_DM = (4*nTimesamples + 3*(nTimesamples>>1))*sizeof(float);
	int nDMs = (free_mem>>1)/memory_required_per_DM;
	printf("DEBUG: nTimesamples: %d; nDMs: %d;\n", nTimesamples, nDMs);
	
	// HOST variables
	int nBoxcars = 32;
	int decimated_timesamples = (nTimesamples>>1);
	int Elements_per_block = PD_NTHREADS*2 - nBoxcars;
	int nBlocks = nTimesamples/Elements_per_block;
    int nRest   = nTimesamples - nBlocks*Elements_per_block;
	if(nRest>0) nBlocks++;
	int iteration = 0;
	int dtm = (nTimesamples>>(iteration+1));
	dtm = dtm - (dtm&1);
	printf("DEBUG: Single pulse detection parameters: nTimesamples=%d; nDMs=%d; nBoxcars=%d; decimated_timesamples=%d; nBlocks=%d; dtm=%d; \n", nTimesamples, nDMs, nBoxcars, decimated_timesamples, nBlocks, dtm);
	
	// Kernel configuration
	dim3 gridSize(nBlocks, nDMs, 1);
    dim3 blockSize(PD_NTHREADS, 1, 1);
	
	// Memory allocation on GPU
	size_t full_size = ((size_t)nDMs)*((size_t)nTimesamples);
	size_t half_size = (full_size/2);
	size_t memory_required = (4*full_size + 3*half_size)*sizeof(float);
	printf("full_size = %zu elements; = %0.2f MB\n", full_size, ((double) (full_size*sizeof(float)))/(1024.0*1024.0));
	printf("Total memory required: %zu bytes = %0.2f MB\n", memory_required, ((double) memory_required)/(1024.0*1024.0));
	
	//if(!device_info.request(memory_required)) {
	//	std::cout << "ERROR:  Could not request memory." << std::endl;
	//	return false;
	//}
	
	int error = 0;
	if ( cudaSuccess != cudaMalloc((void **) &d_input,  sizeof(float)*full_size)) error = -1;
	if ( cudaSuccess != cudaMalloc((void **) &d_decimated,  sizeof(float)*half_size)) error = -1;
	if ( cudaSuccess != cudaMalloc((void **) &d_boxcar_values,  sizeof(float)*2*half_size)) error = -1;
	if ( cudaSuccess != cudaMalloc((void **) &d_MSD_interpolated,  sizeof(float)*2*nBoxcars)) error = -1;
	if ( cudaSuccess != cudaMalloc((void **) &d_output_SNR, sizeof(float)*2*full_size)) error = -1;
	if ( cudaSuccess != cudaMalloc((void **) &d_output_taps, sizeof(ushort)*2*full_size)) error = -1;
	if (error!=0) {
		std::cout << "Error while allocating memory." << std::endl;
		return(0);
	}
	
	// Memory allocation on HOST
	float  *h_input, *h_boxcar_values, *h_decimated, *h_output_SNR, *h_MSD_interpolated;
	ushort *h_output_taps;
	h_input			    = (float *)malloc(full_size*sizeof(float));
	h_boxcar_values     = (float *)malloc(2*half_size*sizeof(float));
	h_decimated		    = (float *)malloc(half_size*sizeof(float));
	h_output_SNR	    = (float *)malloc(2*full_size*sizeof(float));
	h_MSD_interpolated  = (float *)malloc(2*nBoxcars*sizeof(float));
	h_output_taps	    = (ushort *)malloc(2*full_size*sizeof(ushort));
	
	srand(time(NULL));
	for(int d=0; d<nDMs; d++){
		for(int s=0; s<nTimesamples; s++){
			size_t pos = ((size_t) d)*nTimesamples + ((size_t) s);
			//h_input[pos] = rand()/(float)RAND_MAX;
			h_input[pos] = ((int) rand())%10;
		}
	}
	for(int f=0; f<nBoxcars; f++) {
		h_MSD_interpolated[f*2] = 0;
		h_MSD_interpolated[f*2 + 1] = sqrt(f+1);
	}
		
		
	error = 0;
	if ( cudaSuccess != cudaMemcpy(d_input, h_input, full_size*sizeof(float), cudaMemcpyHostToDevice)) error=-1;
	if ( cudaSuccess != cudaMemcpy(d_MSD_interpolated, h_MSD_interpolated, 2*nBoxcars*sizeof(float), cudaMemcpyHostToDevice)) error = -1;
	if(error<0) {
		std::cout << "Error while copying to device." << std::endl;
		return(0);
	}
	
	// GPU single pulse
	call_kernel_SPDT_GPU_1st_plane(gridSize, blockSize, d_input, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, (float2 *) d_MSD_interpolated, nTimesamples, nBoxcars, dtm);
	
	error = 0;
	if ( cudaSuccess != cudaMemcpy(h_boxcar_values, d_boxcar_values, full_size*sizeof(float), cudaMemcpyDeviceToHost)) error=-1;
	if ( cudaSuccess != cudaMemcpy(h_decimated, d_decimated, half_size*sizeof(float), cudaMemcpyDeviceToHost)) error = -1;
	if ( cudaSuccess != cudaMemcpy(h_output_SNR, d_output_SNR, 2*full_size*sizeof(float), cudaMemcpyDeviceToHost)) error = -1;
	if ( cudaSuccess != cudaMemcpy(h_output_taps, d_output_taps, 2*full_size*sizeof(ushort), cudaMemcpyDeviceToHost)) error = -1;
	if(error<0) {
		std::cout << "Error while copying to device." << std::endl;
		return(0);
	}	
	
	int pass = check_boxcar_results(h_input, h_MSD_interpolated, h_boxcar_values, h_decimated, h_output_SNR, h_output_taps, nBoxcars, (size_t) nTimesamples, (size_t) nDMs);
	if(pass==0) std::cout << "Test not passed" << std::endl;
	else std::cout << "Test OK. Passed" << std::endl;

	// deallocation
	free(h_input);
	free(h_boxcar_values);
	free(h_decimated);
	free(h_output_SNR);
	free(h_output_taps);
	
	cudaFree(d_input);
	cudaFree(d_boxcar_values);
	cudaFree(d_decimated);
	cudaFree(d_output_SNR);
	cudaFree(d_output_taps);
	cudaFree(d_MSD_interpolated);
	return 0;
}











