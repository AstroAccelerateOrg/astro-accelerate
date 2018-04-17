#include <helper_cuda.h>

#include "headers/device_MSD_Configuration.h"

#include "device_MSD_shared_kernel_functions.cu"
#include "device_MSD_normal_kernel.cu"
#include "device_MSD_outlier_rejection_kernel.cu"

#include <vector>

//#define MSD_DEBUG

// TODO:
// Remove MSD_legacy


void MSD_init(void) {
	//---------> Specific nVidia stuff
	cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
	cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
}


//---------------------------------------------------------------
//------------- MSD without outlier rejection

int MSD_normal(float *d_MSD, float *d_input, float *d_temp, MSD_Configuration *MSD_conf, cudaStream_t streams) {
	
	#ifdef MSD_DEBUG
	MSD_conf->print();
	#endif

	MSD_init();
	MSD_GPU_limited<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize,0,streams>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset);
	MSD_GPU_final_regular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize,0,streams>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
	
	#ifdef MSD_DEBUG
	float h_MSD[MSD_PARTIAL_SIZE];
	checkCudaErrors(cudaMemcpyAsync(h_MSD, d_MSD, MSD_PARTIAL_SIZE*sizeof(float), cudaMemcpyDeviceToHost, streams)); 
	printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif

	return (0);
}

int MSD_normal(float *d_MSD, float *d_input, int nTimesamples, int nDMs, int offset, cudaStream_t streams){
	int result;
	MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
	float *d_temp;
	checkCudaErrors(cudaMalloc((void **) &d_temp, conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float)));
	result = MSD_normal(d_MSD, d_input, d_temp, &conf, streams);
	checkCudaErrors(cudaFree(d_temp));
	return(result);
}



int MSD_normal_continuous(float *d_MSD, float *d_input, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf, cudaStream_t streams) {

	#ifdef MSD_DEBUG
	MSD_conf->print();
	#endif

	MSD_init();
	MSD_GPU_limited<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize,0,streams>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset);
	MSD_GPU_final_regular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize,0,streams>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, d_previous_partials, MSD_conf->nBlocks_total);
	
	#ifdef MSD_DEBUG
	float h_MSD[3];
	cudaMemcpyAsync(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost, streams); 
	printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif

	return (0);
}

int MSD_normal_continuous(float *d_MSD, float *d_input, float *d_previous_partials, int nTimesamples, int nDMs, int offset, cudaStream_t streams){
	int result;
	MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
	float *d_temp;
	cudaMalloc((void **) &d_temp, conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float));
	result = MSD_normal_continuous(d_MSD, d_input, d_previous_partials, d_temp, &conf, streams);
	cudaFree(d_temp);
	return(result);
}

//------------- MSD without outlier rejection
//---------------------------------------------------------------


//---------------------------------------------------------------
//------------- MSD with outlier rejection

//MSD_BLN_pw
int MSD_outlier_rejection(float *d_MSD, float *d_input, float *d_temp, MSD_Configuration *MSD_conf, float OR_sigma_multiplier, cudaStream_t streams){
	#ifdef MSD_DEBUG
	float h_MSD[3];
	MSD_conf->print();
	#endif
	
	MSD_init();
	MSD_GPU_limited<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize,0,streams>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset);
	MSD_GPU_final_regular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize,0,streams>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
	#ifdef MSD_DEBUG
	cudaMemcpyAsync(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost, streams); 
	printf("Before outlier rejection: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif
	for(int i=0; i<5; i++){
		MSD_BLN_pw_rejection_normal<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize,0,streams>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD,  MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset, OR_sigma_multiplier);
		MSD_GPU_final_nonregular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize,0,streams>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
		#ifdef MSD_DEBUG
		cudaMemcpyAsync(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost, streams); 
		printf("Rejection %d: Mean: %e, Standard deviation: %e; Elements:%zu;\n", i, h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
		printf("---------------------------<\n");
		#endif
	}
	
	#ifdef MSD_DEBUG
	cudaMemcpyAsync(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost, streams); 
	printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif
	
	return(0);
}


int MSD_outlier_rejection(float *d_MSD, float *d_input, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, cudaStream_t streams) {
	int result;
	MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
	float *d_temp;
	cudaMalloc((void **) &d_temp, conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float));
	result = MSD_outlier_rejection(d_MSD, d_input, d_temp, &conf, OR_sigma_multiplier, streams);
	cudaFree(d_temp);
	return(result);
}


//MSD_BLN_pw_continuous
int MSD_outlier_rejection_continuous(float *d_MSD, float *d_input, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf, float OR_sigma_multiplier, cudaStream_t streams){
	#ifdef MSD_DEBUG
	float h_MSD[3];
	MSD_conf->print();
	#endif
	
	MSD_init();
	MSD_GPU_limited<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize,0,streams>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset);
	MSD_GPU_final_regular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize,0,streams>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
	#ifdef MSD_DEBUG
	cudaMemcpyAsync(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost, streams); 
	printf("Before outlier rejection: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif
	for(int i=0; i<5; i++){
		MSD_BLN_pw_rejection_normal<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize,0,streams>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD,  MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset, OR_sigma_multiplier);
		MSD_GPU_final_nonregular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize,0,streams>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
		#ifdef MSD_DEBUG
		cudaMemcpyAsync(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost, streams); 
		printf("Rejection %d: Mean: %e, Standard deviation: %e; Elements:%zu;\n", i, h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
		printf("---------------------------<\n");
		#endif
	}
	MSD_GPU_final_nonregular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize,0,streams>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, d_previous_partials, MSD_conf->nBlocks_total);
	
	#ifdef MSD_DEBUG
	cudaMemcpyAsync(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost, streams); 
	printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif
	
	return(0);
}


int MSD_outlier_rejection_continuous(float *d_MSD, float *d_input, float *d_previous_partials, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, cudaStream_t streams) {
	int result;
	MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
	float *d_temp;
	cudaMalloc((void **) &d_temp, conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float));
	result = MSD_outlier_rejection_continuous(d_MSD, d_input, d_previous_partials, d_temp, &conf, OR_sigma_multiplier, streams);
	cudaFree(d_temp);
	return(result);
}


//MSD_BLN_pw_continuous_OR
int MSD_outlier_rejection_grid(float *d_MSD, float *d_input, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf, float OR_sigma_multiplier, cudaStream_t streams){
	#ifdef MSD_DEBUG
	float h_MSD[3];
	MSD_conf->print();
	#endif
	
	MSD_init();
	MSD_GPU_limited<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize,0,streams>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset);
	MSD_GPU_final_regular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize,0,streams>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
	#ifdef MSD_DEBUG
	cudaMemcpyAsync(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost, streams); 
	printf("Before outlier rejection: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif
	for(int i=0; i<5; i++){
		MSD_BLN_pw_rejection_normal<<<MSD_conf->partials_gridSize,MSD_conf->partials_blockSize,0,streams>>>(d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD,  MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset, OR_sigma_multiplier);
		MSD_GPU_final_nonregular<<<MSD_conf->final_gridSize,MSD_conf->final_blockSize,0,streams>>>(&d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
		#ifdef MSD_DEBUG
		cudaMemcpyAsync(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost, streams); 
		printf("Rejection %d: Mean: %e, Standard deviation: %e; Elements:%zu;\n", i, h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
		printf("---------------------------<\n");
		#endif
	}
	
	#ifdef MSD_DEBUG
	cudaMemcpyAsync(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost, streams); 
	printf("Before grid rejection: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif
	
	MSD_BLN_grid_outlier_rejection_new<<<MSD_conf->final_gridSize, MSD_conf->final_blockSize,0,streams>>>(d_temp, d_MSD, MSD_conf->nBlocks_total+MSD_conf->address, OR_sigma_multiplier);
	
	#ifdef MSD_DEBUG
	cudaMemcpyAsync(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost, streams); 
	printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
	printf("---------------------------<\n");
	#endif
	
	return(0);
}

//------------- MSD with outlier rejection
//---------------------------------------------------------------




void Find_MSD(float *d_MSD, float *d_input, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, int enable_outlier_rejection, cudaStream_t streams){
	MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
	float *d_MSD_workarea;
	checkCudaErrors(cudaMalloc((void **) &d_MSD_workarea, conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float)));
	if(enable_outlier_rejection){
		MSD_outlier_rejection(d_MSD, d_input, d_MSD_workarea, &conf, OR_sigma_multiplier, streams);
	}
	else {
		MSD_normal(d_MSD, d_input, d_MSD_workarea, &conf, streams);
	}
	checkCudaErrors(cudaFree(d_MSD_workarea));
}

void Find_MSD(float *d_MSD, float *d_input, float *d_MSD_workarea, MSD_Configuration *conf, float OR_sigma_multiplier, int enable_outlier_rejection, cudaStream_t streams){
	if(enable_outlier_rejection){
		MSD_outlier_rejection(d_MSD, d_input, d_MSD_workarea, conf, OR_sigma_multiplier, streams);
	}
	else {
		MSD_normal(d_MSD, d_input, d_MSD_workarea, conf, streams);
	}
}

void Find_MSD_continuous(float *d_MSD, float *d_input, float *d_previous_partials, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, int enable_outlier_rejection, cudaStream_t streams){
	MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
	float *d_MSD_workarea;
	checkCudaErrors(cudaMalloc((void **) &d_MSD_workarea, conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float)));
	if(enable_outlier_rejection){
		MSD_outlier_rejection_continuous(d_MSD, d_input, d_previous_partials, d_MSD_workarea, &conf, OR_sigma_multiplier, streams);
	}
	else {
		MSD_normal_continuous(d_MSD, d_input, d_previous_partials, d_MSD_workarea, &conf, streams);
	}
	checkCudaErrors(cudaFree(d_MSD_workarea));
}

void Find_MSD_continuous(float *d_MSD, float *d_input, float *d_previous_partials, float *d_MSD_workarea, MSD_Configuration *conf, float OR_sigma_multiplier, int enable_outlier_rejection, cudaStream_t streams){
	if(enable_outlier_rejection){
		MSD_outlier_rejection_continuous(d_MSD, d_input, d_previous_partials, d_MSD_workarea, conf, OR_sigma_multiplier, streams);
	}
	else {
		MSD_normal_continuous(d_MSD, d_input, d_previous_partials, d_MSD_workarea, conf, streams);
	}
}



//---------------------------------------------------------------
//------------- MSD with outlier rejection on grid

int MSD_grid_outlier_rejection(float *d_MSD, float *d_input, int CellDim_x, int CellDim_y, int nTimesamples, int nDMs, int offset, float multiplier, cudaStream_t streams){
	//---------> Task specific
	int GridSize_x, GridSize_y, x_steps, y_steps, nThreads;
	GridSize_x=(nTimesamples-offset)/CellDim_x;
	GridSize_y=nDMs/CellDim_y;
	x_steps=CellDim_x/WARP;
	if(CellDim_y<HALF_WARP) {
		y_steps  = 1;
		nThreads = WARP*CellDim_y;
	}
	else {
		nThreads = WARP*HALF_WARP;
		y_steps  = CellDim_y/HALF_WARP;
	}

	//---------> Initial phase
	dim3 gridSize(GridSize_x, GridSize_y, 1);
	dim3 blockSize(nThreads, 1, 1);

	//---------> Final phase
	dim3 final_gridSize(1, 1, 1);
	dim3 final_blockSize(WARP*WARP, 1, 1);

	//---------> Allocation of temporary memory
	float *d_output;
	cudaMalloc((void **) &d_output, GridSize_x*GridSize_y*3*sizeof(float));

	//---------> MSD
	MSD_init();
	MSD_BLN_grid_calculate_partials<<<gridSize,blockSize,nThreads*8,streams>>>(d_input, d_output, x_steps, y_steps, nTimesamples, 0);
	MSD_BLN_grid_outlier_rejection<<<final_gridSize, final_blockSize,0,streams>>>(d_output, d_MSD, GridSize_x*GridSize_y, (float) (CellDim_x*CellDim_y), multiplier);

	//---------> De-allocation of temporary memory
	cudaFree(d_output);
	
	return(1);
}

//------------- MSD with outlier rejection on grid
//---------------------------------------------------------------

