//Added by Karel Adamek
//#define SPS_LONG_DEBUG

#include <vector>

#include "headers/params.h"
#include "headers/device_BC_plan.h"
#include "headers/device_MSD_BLN_pw.h"
//#include "headers/device_MSD_BLN_pw_dp.h"
#include "headers/device_MSD_limited.h"
#include "device_SPS_long_kernel.cu"


size_t Get_memory_requirement_of_SPS(){
	return((size_t) (5.5*sizeof(float) + 2*sizeof(ushort)));
}

void Assign_parameters(int f, std::vector<PulseDetection_plan> *PD_plan, int *decimated_timesamples, int *dtm, int *iteration, int *nBoxcars, int *nBlocks, int *output_shift, int *shift, int *startTaps, int *unprocessed_samples, int *total_ut){
	*decimated_timesamples = PD_plan->operator[](f).decimated_timesamples;
	*dtm                   = PD_plan->operator[](f).dtm;
	*iteration             = PD_plan->operator[](f).iteration;
	*nBoxcars              = PD_plan->operator[](f).nBoxcars;
	*nBlocks               = PD_plan->operator[](f).nBlocks;
	*output_shift          = PD_plan->operator[](f).output_shift;
	*shift                 = PD_plan->operator[](f).shift;           
	*startTaps             = PD_plan->operator[](f).startTaps; 
	*unprocessed_samples   = PD_plan->operator[](f).unprocessed_samples;
	*total_ut              = PD_plan->operator[](f).total_ut;	
}

void PD_SEARCH_LONG_init() {
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeEightByte);
}


int PD_SEARCH_LONG_BLN(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nDMs, int nTimesamples) {
	//---------> Task specific
	
	//---------> CUDA block and CUDA grid parameters
	dim3 gridSize(1, 1, 1);
	dim3 blockSize(PD_NTHREADS, 1, 1);
	
	//---------> Pulse detection FIR
	PD_SEARCH_LONG_init();
	
	int f;
	int decimated_timesamples, dtm, iteration, nBoxcars, nBlocks, output_shift, shift, startTaps, unprocessed_samples, total_ut;
	
	// ----------> First iteration
	Assign_parameters(0, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
	gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
	blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
	#ifdef SPS_LONG_DEBUG
	printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples, dtm, iteration ,nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
	#endif
	if(nBlocks>0) PD_GPU_1st_BLN<<<gridSize,blockSize>>>( d_input, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, decimated_timesamples, nBoxcars, dtm);
	
	
	for(f=1; f<max_iteration; f++){
		Assign_parameters(f, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
		gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
		blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
		#ifdef SPS_LONG_DEBUG
		printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples, dtm, iteration, nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
		#endif
		if( (f%2) == 0 ) {
			if(nBlocks>0) PD_GPU_Nth_BLN<<<gridSize,blockSize>>>(&d_input[shift], &d_boxcar_values[nDMs*(nTimesamples>>1)], d_boxcar_values, d_decimated, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
		}
		else {
			if(nBlocks>0) PD_GPU_Nth_BLN<<<gridSize,blockSize>>>(&d_decimated[shift], d_boxcar_values, &d_boxcar_values[nDMs*(nTimesamples>>1)], d_input, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
		}
	}

	return(0);
}


// TODO: this also needs modified kernel since number of tams must be kept local (taps) and not global (star_taps+taps)
int PD_SEARCH_LONG_BLN_EACH(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nDMs, int nTimesamples, float sigma_constant) {
	//---------> Task specific
	float *d_MSD_BV, *d_MSD_DIT;
	if ( cudaSuccess != cudaMalloc((void **) &d_MSD_BV, sizeof(float)*3))  {printf("Allocation error!\n"); exit(1001);}
	if ( cudaSuccess != cudaMalloc((void **) &d_MSD_DIT, sizeof(float)*3)) {printf("Allocation error!\n"); exit(1001);}
	
	#ifdef SPS_LONG_DEBUG
	float h_MSD_BV[3], h_MSD_DIT[3];
	#endif
	
	//---------> CUDA block and CUDA grid parameters
	dim3 gridSize(1, 1, 1);
	dim3 blockSize(PD_NTHREADS, 1, 1);
	
	//---------> Pulse detection FIR
	PD_SEARCH_LONG_init();
	
	int f;
	int decimated_timesamples, dtm, iteration, nBoxcars, nBlocks, output_shift, shift, startTaps, unprocessed_samples, total_ut;
	
	// ----------> First iteration
	Assign_parameters(0, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
	gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
	blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
	//Note: Musim udelat dve SNR jeden pro BV_in and dalsi pro decimated values. Celkove rms je pak rms(BV) + sqrt(ntaps)*rms(decimated)
	MSD_BLN_pw(d_input, d_MSD_BV, nDMs, decimated_timesamples, 0, sigma_constant);
	#ifdef SPS_LONG_DEBUG
	cudaMemcpy(h_MSD_BV, d_MSD_BV, 3*sizeof(float), cudaMemcpyDeviceToHost);
	printf("   decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples, dtm, iteration ,nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
	printf("   MSD BLN point-wise: Mean: %f, Stddev: %f, modifier: %f\n", h_MSD_BV[0], h_MSD_BV[1], h_MSD_BV[2]);
	#endif
	
	if(nBlocks>0) PD_GPU_1st_BLN<<<gridSize,blockSize>>>( d_input, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD_BV, decimated_timesamples, nBoxcars, dtm);
	
	
	for(f=1; f<max_iteration; f++){
		Assign_parameters(f, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
		gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
		blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
		#ifdef SPS_LONG_DEBUG
		printf("   decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples, dtm, iteration, nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
		#endif
		if( (f%2) == 0 ) {
			MSD_BLN_pw(d_input, d_MSD_DIT, nDMs, decimated_timesamples, 0, sigma_constant);
			MSD_BLN_pw(&d_boxcar_values[nDMs*(nTimesamples>>1)], d_MSD_BV, nDMs, decimated_timesamples, PD_plan->operator[](f-1).unprocessed_samples, sigma_constant);
			if(nBlocks>0) PD_GPU_Nth_BLN_EACH<<<gridSize,blockSize>>>(&d_input[shift], &d_boxcar_values[nDMs*(nTimesamples>>1)], d_boxcar_values, d_decimated, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD_BV, d_MSD_DIT, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
		}
		else {
			MSD_BLN_pw(d_decimated, d_MSD_DIT, nDMs, decimated_timesamples, 0, sigma_constant);
			MSD_BLN_pw(d_boxcar_values, d_MSD_BV, nDMs, decimated_timesamples, PD_plan->operator[](f-1).unprocessed_samples, sigma_constant);
			if(nBlocks>0) PD_GPU_Nth_BLN_EACH<<<gridSize,blockSize>>>(&d_decimated[shift], d_boxcar_values, &d_boxcar_values[nDMs*(nTimesamples>>1)], d_input, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD_BV, d_MSD_DIT, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
		}
		
		#ifdef SPS_LONG_DEBUG
		cudaMemcpy(h_MSD_BV, d_MSD_BV, 3*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(h_MSD_DIT, d_MSD_DIT, 3*sizeof(float), cudaMemcpyDeviceToHost);
		printf("   MSD BV point-wise: Mean: %f, Stddev: %f, modifier: %f\n", h_MSD_BV[0], h_MSD_BV[1], h_MSD_BV[2]);
		printf("   MSD DIT point-wise: Mean: %f, Stddev: %f, modifier: %f\n", h_MSD_DIT[0], h_MSD_DIT[1], h_MSD_DIT[2]);
		#endif
	}

	
	cudaFree(d_MSD_BV);
	cudaFree(d_MSD_DIT);
	return(0);
}





int PD_SEARCH_LONG_LINAPPROX(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, float *d_MSD, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nDMs, int nTimesamples) {
	//---------> Task specific
	
	//---------> CUDA block and CUDA grid parameters
	dim3 gridSize(1, 1, 1);
	dim3 blockSize(PD_NTHREADS, 1, 1);
	
	//---------> Pulse detection FIR
	PD_SEARCH_LONG_init();
	
	int f;
	int decimated_timesamples, dtm, iteration, nBoxcars, nBlocks, output_shift, shift, startTaps, unprocessed_samples, total_ut;
	
	// ----------> First iteration
	Assign_parameters(0, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
	gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
	blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
	#ifdef SPS_LONG_DEBUG
	printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples, dtm, iteration ,nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
	#endif
	if(nBlocks>0) PD_GPU_1st_LA<<<gridSize,blockSize>>>( d_input, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, decimated_timesamples, nBoxcars, dtm);
	
	
	for(f=1; f<max_iteration; f++){
		Assign_parameters(f, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
		gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
		blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
		#ifdef SPS_LONG_DEBUG
		printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples, dtm, iteration, nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
		#endif
		if( (f%2) == 0 ) {
			if(nBlocks>0) PD_GPU_Nth_LA<<<gridSize,blockSize>>>(&d_input[shift], &d_boxcar_values[nDMs*(nTimesamples>>1)], d_boxcar_values, d_decimated, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
		}
		else {
			if(nBlocks>0) PD_GPU_Nth_LA<<<gridSize,blockSize>>>(&d_decimated[shift], d_boxcar_values, &d_boxcar_values[nDMs*(nTimesamples>>1)], d_input, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
		}
	}

	return(0);
}


// TODO: this also needs modified kernel since number of tams must be kept local (taps) and not global (star_taps+taps)
int PD_SEARCH_LONG_LINAPPROX_EACH(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nDMs, int nTimesamples) {
	//---------> Task specific
	float *d_MSD, *d_MSD_Nth;
	if ( cudaSuccess != cudaMalloc((void **) &d_MSD, sizeof(float)*4))  {printf("Allocation error!\n"); exit(1001);}
	if ( cudaSuccess != cudaMalloc((void **) &d_MSD_Nth, sizeof(float)*4))  {printf("Allocation error!\n"); exit(1001);}
	#ifdef SPS_LONG_DEBUG
	float h_MSD[4];
	#endif
	
	//---------> CUDA block and CUDA grid parameters
	dim3 gridSize(1, 1, 1);
	dim3 blockSize(PD_NTHREADS, 1, 1);
	
	//---------> Pulse detection FIR
	PD_SEARCH_LONG_init();
	
	int f;
	int decimated_timesamples, dtm, iteration, nBoxcars, nBlocks, output_shift, shift, startTaps, unprocessed_samples, total_ut, error;
	
	// ----------> First iteration
	Assign_parameters(0, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
	gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
	blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
	MSD_linear_approximation(d_input, d_MSD, nBoxcars, nDMs, decimated_timesamples, 0);
	#ifdef SPS_LONG_DEBUG
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost);
	printf("     MSD linear approximation: Mean: %f, Stddev: %f, modifier: %f\n", h_MSD[0], h_MSD[1], h_MSD[2]);
	printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples, dtm, iteration ,nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
	#endif
	if(nBlocks>0) PD_GPU_1st_LA<<<gridSize,blockSize>>>( d_input, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, decimated_timesamples, nBoxcars, dtm);
	
	
	for(f=1; f<max_iteration; f++){
		Assign_parameters(f, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
		gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
		blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
		#ifdef SPS_LONG_DEBUG
		printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples, dtm, iteration, nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
		#endif
		if( (f%2) == 0 ) {
			error = MSD_LA_Nth(&d_input[shift], &d_boxcar_values[nDMs*(nTimesamples>>1)], d_MSD_Nth, d_MSD, nBoxcars, nDMs, decimated_timesamples, 3*unprocessed_samples, (1<<iteration));
			if(nBlocks>0 && error==0) PD_GPU_Nth_LA_EACH<<<gridSize,blockSize>>>(&d_input[shift], &d_boxcar_values[nDMs*(nTimesamples>>1)], d_boxcar_values, d_decimated, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD_Nth, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
			else printf("WARNING: Pulse search does not have enough data!\n");
		}
		else {
			error = MSD_LA_Nth(&d_decimated[shift], d_boxcar_values, d_MSD_Nth, d_MSD, nBoxcars, nDMs, decimated_timesamples, 3*unprocessed_samples, (1<<iteration));
			if(nBlocks>0 && error==0) PD_GPU_Nth_LA_EACH<<<gridSize,blockSize>>>(&d_decimated[shift], d_boxcar_values, &d_boxcar_values[nDMs*(nTimesamples>>1)], d_input, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD_Nth, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
			else printf("WARNING: Pulse search does not have enough data!\n");
		}
		
		#ifdef SPS_LONG_DEBUG
		cudaMemcpy(h_MSD, d_MSD_Nth, 4*sizeof(float), cudaMemcpyDeviceToHost);
		printf("     MSD linear approximation: BV Mean: %f, Stddev: %f, modifier: %f; DIT Mean:%f;\n", h_MSD[0], h_MSD[1], h_MSD[2], h_MSD[3]);
		#endif
	}

	cudaFree(d_MSD);
	return(0);
}




int PD_SEARCH_LONG_BLN_LINAPPROX_EACH(float *d_input, float *d_boxcar_values, float *d_decimated, float *d_output_SNR, ushort *d_output_taps, std::vector<PulseDetection_plan> *PD_plan, int max_iteration, int nDMs, int nTimesamples, float sigma_constant) {
	//---------> Task specific
	float *d_MSD, *d_MSD_Nth;
	if ( cudaSuccess != cudaMalloc((void **) &d_MSD, sizeof(float)*4))  {printf("Allocation error!\n"); exit(1001);}
	if ( cudaSuccess != cudaMalloc((void **) &d_MSD_Nth, sizeof(float)*4))  {printf("Allocation error!\n"); exit(1001);}
	#ifdef SPS_LONG_DEBUG
	float h_MSD[4];
	#endif
	
	//---------> CUDA block and CUDA grid parameters
	dim3 gridSize(1, 1, 1);
	dim3 blockSize(PD_NTHREADS, 1, 1);
	
	//---------> Pulse detection FIR
	PD_SEARCH_LONG_init();
	
	int f;
	int decimated_timesamples, dtm, iteration, nBoxcars, nBlocks, output_shift, shift, startTaps, unprocessed_samples, total_ut;
	
	// ----------> First iteration
	Assign_parameters(0, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
	gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
	blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
	MSD_BLN_LA_pw_normal(d_input, d_MSD, nBoxcars, nDMs, decimated_timesamples, 0, sigma_constant);
	#ifdef SPS_LONG_DEBUG
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost);
	printf("     MSD linear approximation: Mean: %f, Stddev: %f, modifier: %f\n", h_MSD[0], h_MSD[1], h_MSD[2]);
	printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples, dtm, iteration ,nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
	#endif
	if(nBlocks>0) PD_GPU_1st_LA<<<gridSize,blockSize>>>( d_input, d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, decimated_timesamples, nBoxcars, dtm);
	
	
	for(f=1; f<max_iteration; f++){
		Assign_parameters(f, PD_plan, &decimated_timesamples, &dtm, &iteration, &nBoxcars, &nBlocks, &output_shift, &shift, &startTaps, &unprocessed_samples, &total_ut);
		gridSize.x=nBlocks; gridSize.y=nDMs; gridSize.z=1;
		blockSize.x=PD_NTHREADS; blockSize.y=1; blockSize.z=1;
		#ifdef SPS_LONG_DEBUG
		printf("decimated_timesamples:%d; dtm:%d; iteration:%d; nBoxcars:%d; nBlocks:%d; output_shift:%d; shift:%d; startTaps:%d; unprocessed_samples:%d; total_ut:%d;\n",decimated_timesamples, dtm, iteration, nBoxcars ,nBlocks ,output_shift ,shift ,startTaps ,unprocessed_samples ,total_ut);
		#endif
		if( (f%2) == 0 ) {
			MSD_BLN_LA_Nth_pw_normal(&d_input[shift], &d_boxcar_values[nDMs*(nTimesamples>>1)], d_MSD_Nth, d_MSD, nBoxcars, nDMs, decimated_timesamples, 3*unprocessed_samples, (1<<iteration), sigma_constant);
			if(nBlocks>0) PD_GPU_Nth_LA_EACH<<<gridSize,blockSize>>>(&d_input[shift], &d_boxcar_values[nDMs*(nTimesamples>>1)], d_boxcar_values, d_decimated, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD_Nth, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
		}
		else {
			MSD_BLN_LA_Nth_pw_normal(&d_decimated[shift], d_boxcar_values, d_MSD_Nth, d_MSD, nBoxcars, nDMs, decimated_timesamples, 3*unprocessed_samples, (1<<iteration), sigma_constant);
			if(nBlocks>0) PD_GPU_Nth_LA_EACH<<<gridSize,blockSize>>>(&d_decimated[shift], d_boxcar_values, &d_boxcar_values[nDMs*(nTimesamples>>1)], d_input, &d_output_SNR[nDMs*output_shift], &d_output_taps[nDMs*output_shift], d_MSD_Nth, decimated_timesamples, nBoxcars, startTaps, (1<<iteration), dtm);
		}
		
		#ifdef SPS_LONG_DEBUG
		cudaMemcpy(h_MSD, d_MSD_Nth, 4*sizeof(float), cudaMemcpyDeviceToHost);
		printf("     MSD linear approximation: BV Mean: %f, Stddev: %f, modifier: %f; DIT Mean:%f;\n", h_MSD[0], h_MSD[1], h_MSD[2], h_MSD[3]);
		#endif
	}

	cudaFree(d_MSD);
	return(0);
}

