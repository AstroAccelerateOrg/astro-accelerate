//Added by Karel Adamek
//#define MSD_BLN_DEBUG
#define MSD_BLN_END

#include "headers/params.h"
#include "device_MSD_BLN_pw_kernel.cu"

int Choose_divider(int number, int max_divider){
	int seive[12]={2, 3, 4, 5, 7, 11, 13, 17, 19, 23, 29, 31};
	int f, nRest, nBlocks, N, N_accepted;
	
	N=1;N_accepted=1;
	do {
		N=1;
		for(f=0; f<12; f++){
			nBlocks=number/seive[f];
			nRest=number - nBlocks*seive[f];
			if(nRest==0) {
				N=seive[f];
				N_accepted=N_accepted*N;
				break;
			}
		}
		number=number/N;
	} while ( (N_accepted)<=max_divider && N>1 );
	
	return(N_accepted/N);
}

void MSD_BLN_pw_init(){
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
}

int MSD_BLN_pw(float *d_input, float *d_MSD, int nDMs, int nTimesamples, int offset, float bln_sigma_constant){
	//---------> Task specific
	ushort nBlocks_x, nBlocks_y;
	int	nBlocks_total, nSteps_x, nSteps_y, nRest;
	float *d_output;
	
	nBlocks_x=0; nRest=0;
	nSteps_x=MSD_PW_NTHREADS;
	nBlocks_x = nBlocks_x + (nTimesamples-offset)/(nSteps_x);
	nRest=nTimesamples - offset - nBlocks_x*nSteps_x;
	if(nRest>32) nBlocks_x++;
	
	nSteps_y = Choose_divider(nDMs,64);
	nBlocks_y=nDMs/nSteps_y;
	nBlocks_total=nBlocks_x*nBlocks_y;

	int nThreads=2048;
	int itemp=0;
	while(itemp==0 && nThreads>32){
		nThreads=(nThreads>>1);
		itemp=(int) (nBlocks_total/(nThreads*32));
	}
	if(nThreads<32) nThreads=32;
	
	dim3 gridSize(nBlocks_x, nBlocks_y, 1);
	dim3 blockSize(MSD_PW_NTHREADS, 1, 1);
	dim3 final_gridSize(1, 1, 1);
	dim3 final_blockSize(nThreads, 1, 1);
	
	//---------> Allocation of temporary memory
	cudaMalloc((void **) &d_output, nBlocks_total*3*sizeof(float));
	
	MSD_BLN_pw_init();
	MSD_BLN_pw_no_rejection<<<gridSize,blockSize>>>(d_input, d_output, nDMs/nBlocks_y, nTimesamples, offset);
	MSD_GPU_final_regular<<<final_gridSize,final_blockSize>>>(d_output, d_MSD, nBlocks_total);
	for(int i=0; i<5; i++){
		MSD_BLN_pw_rejection_normal<<<gridSize,blockSize>>>(d_input, d_output, d_MSD, nDMs/nBlocks_y, nTimesamples, offset, bln_sigma_constant);
		MSD_GPU_final_nonregular<<<final_gridSize, final_blockSize>>>(d_output, d_MSD, nBlocks_total);
	}
	
	cudaFree(d_output);
	
	return(0);
}

int MSD_BLN_pw(float *d_input, float *d_MSD, float *d_temp, MSD_Configuration *MSD_conf, float bln_sigma_constant){
	
	#ifdef MSD_BLN_DEBUG
	float h_MSD[3];
	MSD_conf->print();
	#endif
	
	MSD_BLN_pw_init();
	MSD_BLN_pw_no_rejection<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset);
	MSD_GPU_final_regular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
	#ifdef MSD_BLN_DEBUG
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	printf("Before outlier rejection: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif
	for(int i=0; i<5; i++){
		MSD_BLN_pw_rejection_normal<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD,  MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset, bln_sigma_constant);
		MSD_GPU_final_nonregular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
		#ifdef MSD_BLN_DEBUG
		cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
		printf("Rejection %d: Mean: %e, Standard deviation: %e; Elements:%zu;\n", i, h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
		printf("---------------------------<\n");
		#endif
	}
	
	#ifdef MSD_BLN_DEBUG
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif
	
	return(0);
}

int MSD_BLN_pw_continuous(float *d_input, float *d_MSD, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf, float bln_sigma_constant){
	
	#ifdef MSD_BLN_DEBUG
	float h_MSD[3];
	MSD_conf->print();
	#endif
	
	MSD_BLN_pw_init();
	MSD_BLN_pw_no_rejection<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset);
	MSD_GPU_final_regular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
	#ifdef MSD_BLN_DEBUG
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	printf("Before outlier rejection: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif
	for(int i=0; i<5; i++){
		MSD_BLN_pw_rejection_normal<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD,  MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset, bln_sigma_constant);
		MSD_GPU_final_nonregular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
		#ifdef MSD_BLN_DEBUG
		cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
		printf("Rejection %d: Mean: %e, Standard deviation: %e; Elements:%zu;\n", i, h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
		printf("---------------------------<\n");
		#endif
	}
	
	MSD_GPU_final_nonregular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, d_previous_partials, MSD_conf->nBlocks_total);
	
	#ifdef MSD_BLN_DEBUG
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	printf("Before grid rejection: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif
	
	return(0);
}

int MSD_BLN_pw_continuous_OR(float *d_input, float *d_MSD, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf, float bln_sigma_constant){
	
	#ifdef MSD_BLN_DEBUG
	float h_MSD[3];
	MSD_conf->print();
	#endif
	
	MSD_BLN_pw_init();
	MSD_BLN_pw_no_rejection<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset);
	MSD_GPU_final_regular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
	#ifdef MSD_BLN_DEBUG
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	printf("Before outlier rejection: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif
	for(int i=0; i<5; i++){
		MSD_BLN_pw_rejection_normal<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD,  MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset, bln_sigma_constant);
		MSD_GPU_final_nonregular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
		#ifdef MSD_BLN_DEBUG
		cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
		printf("Rejection %d: Mean: %e, Standard deviation: %e; Elements:%zu;\n", i, h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
		printf("---------------------------<\n");
		#endif
	}
	
	#ifdef MSD_BLN_DEBUG
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	printf("Before grid rejection: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif
	
	MSD_BLN_grid_outlier_rejection_new<<<MSD_conf->final_gridSize, MSD_conf->final_blockSize>>>(d_temp, d_MSD, MSD_conf->nBlocks_total+MSD_conf->address, bln_sigma_constant);
	
	#ifdef MSD_BLN_DEBUG
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif
	
	return(0);
}

int MSD_BLN_LA_pw_normal(float *d_input, float *d_MSD_T, int nTaps, int nDMs, int nTimesamples, int offset, float bln_sigma_constant){
	//---------> Task specific
	int nBlocks_x, nBlocks_x_LA, nBlocks_y;
	int	nBlocks_total_LA, nBlocks_total, nSteps_x_LA, nSteps_x, nSteps_y, nRest, nRest_LA;
	float *d_output;
	float *d_output_taps;
	float *d_MSD_T_base;
	#ifdef MSD_BLN_END
	float h_MSD_break[3], h_MSD_break_old[3];
	#endif
	
	//----------------> y-step is same for both
	nSteps_y  = Choose_divider(nDMs,64);
	nBlocks_y = nDMs/nSteps_y;
	
	//----------------> LA_kernel
	nSteps_x_LA  = 2*MSD_PW_NTHREADS-nTaps+4;
	nBlocks_x_LA = (int) ((nTimesamples-offset)/nSteps_x_LA);
	nRest_LA     = nTimesamples - offset - nBlocks_x_LA*nSteps_x_LA;
	if(nRest_LA>128) nBlocks_x_LA++;
	nBlocks_total_LA = nBlocks_x_LA*nBlocks_y;
	dim3 gridSize_LA(nBlocks_x_LA, nBlocks_y, 1);
	
	//----------------> Base_kernel
	nSteps_x  = MSD_PW_NTHREADS;
	nBlocks_x = (int) ((nTimesamples-offset)/nSteps_x);
	nRest     = nTimesamples - offset - nBlocks_x*nSteps_x;
	if(nRest>128) nBlocks_x++;
	nBlocks_total = nBlocks_x*nBlocks_y;
	dim3 gridSize(nBlocks_x, nBlocks_y, 1);		
	
	dim3 blockSize(MSD_PW_NTHREADS, 1, 1);
	
	int nThreads=2048;
	int itemp=0;
	while(itemp==0 && nThreads>32){
		nThreads=(nThreads>>1);
		itemp=(int) (nBlocks_total/(nThreads*32));
	}
	if(nThreads<32) nThreads=32;
	
	dim3 final_gridSize(1, 1, 1);
	dim3 final_blockSize(nThreads, 1, 1);
	
	#ifdef MSD_BLN_DEBUG
	float h_MSD_T[4], h_MSD_T_base[3];
	printf("\n\n");
	printf("----------------> MSD debug:\n");
	printf("nDMs:%d; nTimesamples:%d; nTaps:%d; offset:%d\n", nDMs, nTimesamples, nTaps, offset);
	printf("Kernel for calculating partials: (MSD_BLN_LA_pw_normal)\n");
	printf("d_MSD_T kernel\n");
	printf("ThreadBlocks (TB) in x:%d; Elements processed by TB in x:%d; Remainder in x:%d", nBlocks_x_LA, nSteps_x_LA, nRest_LA);
	if(nRest_LA>3*nTaps) printf(" is processed\n");
	else printf(" is not processed\n");
	printf("ThreadBlocks (TB) in y:%d; Elements processed by TB in y:%d; Remainder in y:%d is processed\n", nBlocks_y, nSteps_y, 0);
	printf("d_MSD_T_base kernel\n");
	printf("ThreadBlocks (TB) in x:%d; Elements processed by TB in x:%d; Remainder in x:%d", nBlocks_x, nSteps_x, nRest);
	if(nRest>3*nTaps) printf(" is processed\n");
	else printf(" is not processed\n");
	printf("ThreadBlocks (TB) in y:%d; Elements processed by TB in y:%d; Remainder in y:%d is processed\n", nBlocks_y, nSteps_y, 0);
	printf("gridSize_LA=(%d,%d,%d)\n", gridSize_LA.x, gridSize_LA.y, gridSize_LA.z);
	printf("gridSize=(%d,%d,%d)\n", gridSize.x, gridSize.y, gridSize.z);
	printf("blockSize=(%d,%d,%d)\n", blockSize.x, blockSize.y, blockSize.z);
	printf("Shared memory required: %0.3f B\n", (float) (MSD_PW_NTHREADS*24));
	printf("Kernel for final calculation of mean and standard deviation:\n");
	size_t free_mem,total_mem;
	cudaMemGetInfo(&free_mem,&total_mem);
	printf("Memory required for temporary storage:%0.3f MB which is %d floats\n",(nBlocks_total*3*2*sizeof(float))/(1024.0*1024),  nBlocks_total*3*2);
	printf("Memory available:%0.3f MB \n", ((float) free_mem)/(1024.0*1024.0) );
	printf("gridSize=(%d,%d,%d)\n", final_gridSize.x, final_gridSize.y, final_gridSize.z);
	printf("blockSize=(%d,%d,%d)\n", final_blockSize.x, final_blockSize.y, final_blockSize.z);	
	printf("\n");
	#endif
	
	
	if(nBlocks_total>0 && nBlocks_total_LA>0){
		//---------> Allocation of temporary memory
		cudaMalloc((void **) &d_output, nBlocks_total*3*sizeof(float));
		cudaMalloc((void **) &d_output_taps, nBlocks_total*3*sizeof(float));
		cudaMalloc((void **) &d_MSD_T_base, 3*sizeof(float));
		
		MSD_BLN_pw_init();
		MSD_GPU_LA_ALL_no_rejection<<<gridSize_LA,blockSize>>>(d_input, d_output, d_output_taps, nSteps_y, nTaps, nTimesamples, offset);
		MSD_BLN_pw_no_rejection<<<gridSize,blockSize>>>(d_input, d_output, nSteps_y, nTimesamples, offset);
		MSD_GPU_final_regular<<<final_gridSize, final_blockSize>>>(d_output, d_MSD_T_base, nBlocks_total);
		MSD_GPU_final_regular<<<final_gridSize, final_blockSize>>>(d_output_taps, d_MSD_T, nBlocks_total_LA);
		//MSD_GPU_limited_final_create_linear_approx_regular<<<final_gridSize, final_blockSize>>>(d_output_taps, d_MSD_T, d_MSD_T_base, nTaps, nBlocks_total);
		
		#ifdef MSD_BLN_DEBUG
		cudaMemcpy(h_MSD_T, d_MSD_T, 3*sizeof(float), cudaMemcpyDeviceToHost); 
		cudaMemcpy(h_MSD_T_base, d_MSD_T_base, 3*sizeof(float), cudaMemcpyDeviceToHost);
		printf("d_MSD_T: Mean: %f, StDev: %f; Elements:%f; \n", h_MSD_T[0], h_MSD_T[1], h_MSD_T[2]);
		printf("d_MSD_T_base: Mean: %f, StDev: %f; Elements:%f; \n", h_MSD_T_base[0], h_MSD_T_base[1], h_MSD_T_base[2]);
		#endif
		#ifdef MSD_BLN_END
		cudaMemcpy(h_MSD_break, d_MSD_T, 3*sizeof(float), cudaMemcpyDeviceToHost);
		#endif
		
		for(int i=0; i<5; i++){
			MSD_GPU_LA_ALL_pw_rejection<<<gridSize_LA,blockSize>>>(d_input, d_output, d_output_taps, d_MSD_T, d_MSD_T_base, bln_sigma_constant, nSteps_y, nTaps, nTimesamples, offset);
			MSD_BLN_pw_rejection_normal<<<gridSize,blockSize>>>(d_input, d_output, d_MSD_T_base, nSteps_y, nTimesamples, offset, bln_sigma_constant);
			MSD_GPU_final_nonregular<<<final_gridSize, final_blockSize>>>(d_output, d_MSD_T_base, nBlocks_total);
			MSD_GPU_final_nonregular<<<final_gridSize, final_blockSize>>>(d_output_taps, d_MSD_T, nBlocks_total_LA);
			//MSD_GPU_limited_final_create_linear_approx_nonregular<<<final_gridSize, final_blockSize>>>(d_output_taps, d_MSD_T, d_MSD_T_base, nTaps, nBlocks_total);
			#ifdef MSD_BLN_DEBUG
			printf("-------------------\n");
			cudaMemcpy(h_MSD_T, d_MSD_T, 3*sizeof(float), cudaMemcpyDeviceToHost); 
			cudaMemcpy(h_MSD_T_base, d_MSD_T_base, 3*sizeof(float), cudaMemcpyDeviceToHost);
			printf("d_MSD_T: Mean: %f, StDev: %f; Elements:%f; \n", h_MSD_T[0], h_MSD_T[1], h_MSD_T[2]);
			printf("d_MSD_T_base: Mean: %f, StDev: %f; Elements:%f; \n", h_MSD_T_base[0], h_MSD_T_base[1], h_MSD_T_base[2]);
			#endif
			#ifdef MSD_BLN_END
			h_MSD_break_old[0] = h_MSD_break[0]; h_MSD_break_old[1] = h_MSD_break[1]; h_MSD_break_old[2] = h_MSD_break[2];
			cudaMemcpy(h_MSD_break, d_MSD_T, 3*sizeof(float), cudaMemcpyDeviceToHost);
			if( (h_MSD_break_old[1]-h_MSD_break[1])<1e-4 ) {
				break;
			}
			#endif
		}
		//MSD_GPU_LA_ALL_pw_rejection<<<gridSize,blockSize>>>(d_input, d_output, d_output_taps, d_MSD_T, d_MSD_T_base, bln_sigma_constant, nSteps_y, nTaps, nTimesamples, offset);
		//MSD_BLN_pw_rejection_normal<<<gridSize,blockSize>>>(d_input, d_output, d_MSD_T_base, nSteps_y, nTimesamples, offset, bln_sigma_constant);
		//MSD_GPU_final_nonregular<<<final_gridSize, final_blockSize>>>(d_output, d_MSD_T_base, nBlocks_total);
		//MSD_GPU_final_nonregular<<<final_gridSize, final_blockSize>>>(d_output_taps, d_MSD_T, nBlocks_total);
		MSD_GPU_final_create_LA_nonregular<<<final_gridSize, final_blockSize>>>(d_output_taps, d_MSD_T, d_MSD_T_base, nTaps, nBlocks_total_LA);
		
		#ifdef MSD_BLN_DEBUG
		printf("-------------------\n");
		cudaMemcpy(h_MSD_T, d_MSD_T, 3*sizeof(float), cudaMemcpyDeviceToHost); 
		cudaMemcpy(h_MSD_T_base, d_MSD_T_base, 3*sizeof(float), cudaMemcpyDeviceToHost);
		printf("Final output: Mean base: %f, StDev base: %f; Modifier:%f; \n", h_MSD_T[0], h_MSD_T[1], h_MSD_T[2]);
		printf("d_MSD_T_base: Mean: %f, StDev: %f; Elements:%f; \n", h_MSD_T_base[0], h_MSD_T_base[1], h_MSD_T_base[2]);
		printf("-----------------------------<\n");
		#endif	
		
		cudaFree(d_output);
		cudaFree(d_output_taps);
		cudaFree(d_MSD_T_base);
		//-------------------------------
	}
	else {
		printf("Number of time samples is too small! Increase number of samples send to the boxcar filters. (MSD_BLN_LA_pw_normal)\n");
		exit(1002);
	}
	
	return(0);
}


int MSD_BLN_LA_Nth_pw_normal(float *d_input, float *d_bv_in, float *d_MSD_T, float *d_MSD_DIT, int nTaps, int nDMs, int nTimesamples, int offset, int DIT_value, float bln_sigma_constant){
	//---------> Task specific
	int nBlocks_x, nBlocks_x_LA, nBlocks_y;
	int	nBlocks_total_LA, nBlocks_total, nSteps_x_LA, nSteps_x, nSteps_y, nRest, nRest_LA;
	float *d_output;
	float *d_output_taps;
	float *d_MSD_T_base;
	#ifdef MSD_BLN_END
	float h_MSD_break[3], h_MSD_break_old[3];
	#endif
	
	//----------------> y-step is same for both
	nSteps_y  = Choose_divider(nDMs,64);
	nBlocks_y = nDMs/nSteps_y;
	
	//----------------> LA_kernel
	nSteps_x_LA  = 2*MSD_PW_NTHREADS-nTaps+4;
	nBlocks_x_LA = (int) ((nTimesamples-offset)/nSteps_x_LA);
	nRest_LA     = nTimesamples - offset - nBlocks_x_LA*nSteps_x_LA;
	if(nRest_LA>128) nBlocks_x_LA++;
	nBlocks_total_LA = nBlocks_x_LA*nBlocks_y;
	dim3 gridSize_LA(nBlocks_x_LA, nBlocks_y, 1);
	
	//----------------> Base_kernel
	nSteps_x  = MSD_PW_NTHREADS;
	nBlocks_x = (int) ((nTimesamples-offset)/nSteps_x);
	nRest     = nTimesamples - offset - nBlocks_x*nSteps_x;
	if(nRest>128) nBlocks_x++;
	nBlocks_total = nBlocks_x*nBlocks_y;
	dim3 gridSize(nBlocks_x, nBlocks_y, 1);	
	
	dim3 blockSize(MSD_PW_NTHREADS, 1, 1);
	
	int nThreads=2048;
	int itemp=0;
	while(itemp==0 && nThreads>32){
		nThreads=(nThreads>>1);
		itemp=(int) (nBlocks_total/(nThreads*32));
	}
	if(nThreads<32) nThreads=32;
	
	dim3 final_gridSize(1, 1, 1);
	dim3 final_blockSize(nThreads, 1, 1);
	
	
	#ifdef MSD_BLN_DEBUG
	float h_MSD_T[4], h_MSD_T_base[3];
	printf("\n\n");
	printf("----------------> MSD debug:\n");
	printf("nDMs:%d; nTimesamples:%d; nTaps:%d; offset:%d\n", nDMs, nTimesamples, nTaps, offset);
	printf("Kernel for calculating partials: (MSD_BLN_LA_Nth_pw_normal)\n");
	printf("d_MSD_T kernel\n");
	printf("ThreadBlocks (TB) in x:%d; Elements processed by TB in x:%d; Remainder in x:%d", nBlocks_x_LA, nSteps_x_LA, nRest_LA);
	if(nRest_LA>3*nTaps) printf(" is processed\n");
	else printf(" is not processed\n");
	printf("ThreadBlocks (TB) in y:%d; Elements processed by TB in y:%d; Remainder in y:%d is processed\n", nBlocks_y, nSteps_y, 0);
	printf("d_MSD_T_base kernel\n");
	printf("ThreadBlocks (TB) in x:%d; Elements processed by TB in x:%d; Remainder in x:%d", nBlocks_x, nSteps_x, nRest);
	if(nRest>3*nTaps) printf(" is processed\n");
	else printf(" is not processed\n");
	printf("ThreadBlocks (TB) in y:%d; Elements processed by TB in y:%d; Remainder in y:%d is processed\n", nBlocks_y, nSteps_y, 0);
	printf("gridSize_LA=(%d,%d,%d)\n", gridSize_LA.x, gridSize_LA.y, gridSize_LA.z);
	printf("gridSize=(%d,%d,%d)\n", gridSize.x, gridSize.y, gridSize.z);
	printf("blockSize=(%d,%d,%d)\n", blockSize.x, blockSize.y, blockSize.z);
	printf("Shared memory required: %0.3f B\n", (float) (MSD_PW_NTHREADS*24));
	printf("Kernel for final calculation of mean and standard deviation:\n");
	size_t free_mem,total_mem;
	cudaMemGetInfo(&free_mem,&total_mem);
	printf("Memory required for temporary storage:%0.3f MB which is %d floats\n",(nBlocks_total*3*2*sizeof(float))/(1024.0*1024),  nBlocks_total*3*2);
	printf("Memory available:%0.3f MB \n", ((float) free_mem)/(1024.0*1024.0) );
	printf("gridSize=(%d,%d,%d)\n", final_gridSize.x, final_gridSize.y, final_gridSize.z);
	printf("blockSize=(%d,%d,%d)\n", final_blockSize.x, final_blockSize.y, final_blockSize.z);	
	printf("\n");
	#endif
	
	if(nBlocks_total>0 && nBlocks_total_LA>0){
		//---------> Allocation of temporary memory
		cudaMalloc((void **) &d_output, nBlocks_total*3*sizeof(float));
		cudaMalloc((void **) &d_output_taps, nBlocks_total*3*sizeof(float));
		cudaMalloc((void **) &d_MSD_T_base, 3*sizeof(float));
		
		MSD_BLN_pw_init();
		MSD_GPU_LA_ALL_Nth_no_rejection<<<gridSize_LA,blockSize>>>(d_input, d_bv_in, d_output, d_output_taps, nSteps_y, nTaps, nTimesamples, offset);
		MSD_BLN_pw_no_rejection<<<gridSize,blockSize>>>(d_bv_in, d_output, nSteps_y, nTimesamples, offset);
		MSD_GPU_final_regular<<<final_gridSize, final_blockSize>>>(d_output, d_MSD_T_base, nBlocks_total);
		MSD_GPU_final_regular<<<final_gridSize, final_blockSize>>>(d_output_taps, d_MSD_T, nBlocks_total_LA);
		//MSD_GPU_limited_final_create_linear_approx_regular<<<final_gridSize, final_blockSize>>>(d_output_taps, d_MSD_T, d_MSD_T_base, nTaps, nBlocks_total);
		
		#ifdef MSD_BLN_DEBUG
		cudaMemcpy(h_MSD_T, d_MSD_T, 3*sizeof(float), cudaMemcpyDeviceToHost); 
		cudaMemcpy(h_MSD_T_base, d_MSD_T_base, 3*sizeof(float), cudaMemcpyDeviceToHost);
		printf("d_MSD_T (end):        Mean: %f, StDev: %f; Elements:%f; \n", h_MSD_T[0], h_MSD_T[1], h_MSD_T[2]);
		printf("d_MSD_T_base (start): Mean: %f, StDev: %f; Elements:%f; \n", h_MSD_T_base[0], h_MSD_T_base[1], h_MSD_T_base[2]);
		#endif
		#ifdef MSD_BLN_END
		cudaMemcpy(h_MSD_break, d_MSD_T, 3*sizeof(float), cudaMemcpyDeviceToHost);
		#endif
		for(int i=0; i<5; i++){
			MSD_GPU_LA_ALL_Nth_pw_rejection<<<gridSize_LA,blockSize>>>(d_input, d_bv_in, d_output, d_output_taps, d_MSD_T, d_MSD_T_base, bln_sigma_constant, nSteps_y, nTaps, nTimesamples, offset);
			MSD_BLN_pw_rejection_normal<<<gridSize,blockSize>>>(d_bv_in, d_output, d_MSD_T_base, nSteps_y, nTimesamples, offset, bln_sigma_constant);
			MSD_GPU_final_nonregular<<<final_gridSize, final_blockSize>>>(d_output, d_MSD_T_base, nBlocks_total);
			MSD_GPU_final_nonregular<<<final_gridSize, final_blockSize>>>(d_output_taps, d_MSD_T, nBlocks_total_LA);
			//MSD_GPU_limited_final_create_linear_approx_nonregular<<<final_gridSize, final_blockSize>>>(d_output_taps, d_MSD_T, d_MSD_T_base, nTaps, nBlocks_total);
			#ifdef MSD_BLN_DEBUG
			printf("-------------------\n");
			cudaMemcpy(h_MSD_T, d_MSD_T, 3*sizeof(float), cudaMemcpyDeviceToHost); 
			cudaMemcpy(h_MSD_T_base, d_MSD_T_base, 3*sizeof(float), cudaMemcpyDeviceToHost);
			printf("d_MSD_T (end):        Mean: %f, StDev: %f; Elements:%f; \n", h_MSD_T[0], h_MSD_T[1], h_MSD_T[2]);
			printf("d_MSD_T_base (start): Mean: %f, StDev: %f; Elements:%f; \n", h_MSD_T_base[0], h_MSD_T_base[1], h_MSD_T_base[2]);
			#endif
			#ifdef MSD_BLN_END
			h_MSD_break_old[0] = h_MSD_break[0]; h_MSD_break_old[1] = h_MSD_break[1]; h_MSD_break_old[2] = h_MSD_break[2];
			cudaMemcpy(h_MSD_break, d_MSD_T, 3*sizeof(float), cudaMemcpyDeviceToHost);
			if( (h_MSD_break_old[1]-h_MSD_break[1])<1e-4 ) {
				break;
			}
			#endif
		}
		//MSD_GPU_LA_ALL_Nth_pw_rejection<<<gridSize,blockSize>>>(d_input, d_bv_in, d_output, d_output_taps, d_MSD_T, d_MSD_T_base, bln_sigma_constant, nSteps_y, nTaps, nTimesamples, offset);
		//MSD_BLN_pw_rejection_normal<<<gridSize,blockSize>>>(d_bv_in, d_output, d_MSD_T_base, nSteps_y, nTimesamples, offset, bln_sigma_constant);
		//MSD_GPU_final_nonregular<<<final_gridSize, final_blockSize>>>(d_output, d_MSD_T_base, nBlocks_total);
		MSD_GPU_final_create_LA_Nth_nonregular<<<final_gridSize, final_blockSize>>>(d_output_taps, d_MSD_T, d_MSD_T_base, d_MSD_DIT, nTaps, nBlocks_total_LA, DIT_value);
		//MSD_GPU_limited_final_create_linear_approx_nonregular_final_Nth<<<final_gridSize, final_blockSize>>>(d_output_taps, d_MSD_T, d_MSD_T_base, d_MSD_DIT, nTaps, nBlocks_total, DIT_value);
		
		#ifdef MSD_BLN_DEBUG
		printf("-------------------\n");
		cudaMemcpy(h_MSD_T, d_MSD_T, 3*sizeof(float), cudaMemcpyDeviceToHost); 
		cudaMemcpy(h_MSD_T_base, d_MSD_T_base, 3*sizeof(float), cudaMemcpyDeviceToHost);
		printf("Final output: Mean base: %f, StDev base: %f; Modifier:%f; \n", h_MSD_T[0], h_MSD_T[1], h_MSD_T[2]);
		printf("d_MSD_T_base: Mean: %f, StDev: %f; Elements:%f; \n", h_MSD_T_base[0], h_MSD_T_base[1], h_MSD_T_base[2]);
		printf("-----------------------------<\n");
		#endif
		
		cudaFree(d_output);
		cudaFree(d_output_taps);
		cudaFree(d_MSD_T_base);
		//-------------------------------
	}
	else {
		printf("Number of time samples is too small! Increase number of samples send to the boxcar filters. (MSD_BLN_LA_Nth_pw_normal)\n");
		exit(1002);
	}
	
	return(0);
}
