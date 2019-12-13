#include "aa_device_MSD.hpp"

#include <string>
#include "aa_log.hpp"

//#define MSD_DEBUG

// \todo Remove MSD_legacy

namespace astroaccelerate {


void MSD_init(void) {
  //---------> Specific nVidia stuff
  cudaDeviceSetCacheConfig (cudaFuncCachePreferShared);
  cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
}


//---------------------------------------------------------------
//------------- MSD without outlier rejection

int MSD_normal(float *d_MSD, float *d_input, float *d_temp, MSD_Configuration *MSD_conf) {
	
#ifdef MSD_DEBUG
  MSD_conf->print();
#endif

  MSD_init();
  call_kernel_MSD_GPU_limited(MSD_conf->partials_gridSize,MSD_conf->partials_blockSize, d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int)MSD_conf->offset);
  call_kernel_MSD_GPU_final_regular(MSD_conf->final_gridSize,MSD_conf->final_blockSize,
				    &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
	
#ifdef MSD_DEBUG
  float h_MSD[MSD_PARTIAL_SIZE];
  cudaError_t e = cudaMemcpy(h_MSD, d_MSD, MSD_PARTIAL_SIZE*sizeof(float), cudaMemcpyDeviceToHost);

  if(e != cudaSuccess) {
    LOG(log_level::error, "Could not cudaMemcpy in aa_device_MSD.cu (" + std::string(cudaGetErrorString(e)) + ")");
  }
  
  printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
  printf("---------------------------<\n");
#endif

  return (0);
}

int MSD_normal(float *d_MSD, float *d_input, int nTimesamples, int nDMs, int offset){
  int result;
  MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
  float *d_temp;
  cudaError_t e = cudaMalloc((void **) &d_temp, conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float));

  if(e != cudaSuccess) {
    LOG(log_level::error, "Could not cudaMemcpy in aa_device_MSD.cu (" + std::string(cudaGetErrorString(e)) + ")");
  }
  
  result = MSD_normal(d_MSD, d_input, d_temp, &conf);
  e = cudaFree(d_temp);
  
  if(e != cudaSuccess) {
    LOG(log_level::error, "Could not cudaMemcpy in aa_device_MSD.cu (" + std::string(cudaGetErrorString(e)) + ")");
  }
  
  return(result);
}



int MSD_normal_continuous(float *d_MSD, float *d_input, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf) {

#ifdef MSD_DEBUG
  MSD_conf->print();
#endif

  MSD_init();
  call_kernel_MSD_GPU_limited(MSD_conf->partials_gridSize,MSD_conf->partials_blockSize, d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int)MSD_conf->offset);
  call_kernel_MSD_GPU_final_regular(MSD_conf->final_gridSize, MSD_conf->final_blockSize,
				    &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE],
				    d_MSD, d_previous_partials, MSD_conf->nBlocks_total);
	
#ifdef MSD_DEBUG
  float h_MSD[3];
  cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
  printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
  printf("---------------------------<\n");
#endif

  return (0);
}

int MSD_normal_continuous(float *d_MSD, float *d_input, float *d_previous_partials, int nTimesamples, int nDMs, int offset){
  int result;
  MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
  float *d_temp;
  cudaMalloc((void **) &d_temp, conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float));
  result = MSD_normal_continuous(d_MSD, d_input, d_previous_partials, d_temp, &conf);
  cudaFree(d_temp);
  return(result);
}

//------------- MSD without outlier rejection
//---------------------------------------------------------------


//---------------------------------------------------------------
//------------- MSD with outlier rejection

//MSD_BLN_pw
int MSD_outlier_rejection(float *d_MSD, float *d_input, float *d_temp, MSD_Configuration *MSD_conf, float OR_sigma_multiplier){
	#ifdef MSD_DEBUG
		float h_MSD[3];
		MSD_conf->print();
	#endif
	
	MSD_init();
	call_kernel_MSD_GPU_calculate_partials_2d_and_minmax(MSD_conf->partials_gridSize,MSD_conf->partials_blockSize, d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int)MSD_conf->offset);
	call_kernel_MSD_GPU_final_regular(MSD_conf->final_gridSize, MSD_conf->final_blockSize, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);

	#ifdef MSD_DEBUG
		cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
		printf("Before outlier rejection: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
		printf("---------------------------<\n");
	#endif
	for(int i=0; i<5; i++){
		call_kernel_MSD_BLN_calculate_partials_2d_and_minmax_with_outlier_rejection(MSD_conf->partials_gridSize,MSD_conf->partials_blockSize, d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD,  MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset, OR_sigma_multiplier);
		call_kernel_MSD_GPU_final_nonregular(MSD_conf->final_gridSize,MSD_conf->final_blockSize, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD, MSD_conf->nBlocks_total);
		#ifdef MSD_DEBUG
			cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
			printf("Rejection %d: Mean: %e, Standard deviation: %e; Elements:%zu;\n", i, h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
			printf("---------------------------<\n");
		#endif
	}
	
	#ifdef MSD_DEBUG
		cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
		printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
		printf("---------------------------<\n");
	#endif
	
	return(0);
}


int MSD_outlier_rejection(float *d_MSD, float *d_input, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier) {
  int result;
  MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
  float *d_temp;
  cudaMalloc((void **) &d_temp, conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float));
  result = MSD_outlier_rejection(d_MSD, d_input, d_temp, &conf, OR_sigma_multiplier);
  cudaFree(d_temp);
  return(result);
}


//MSD_BLN_pw_continuous
int MSD_outlier_rejection_continuous(float *d_MSD, float *d_input, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf, float OR_sigma_multiplier){
#ifdef MSD_DEBUG
  float h_MSD[3];
  MSD_conf->print();
#endif
	
  MSD_init();
  call_kernel_MSD_GPU_calculate_partials_2d_and_minmax(MSD_conf->partials_gridSize,MSD_conf->partials_blockSize, d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int)MSD_conf->offset);
  call_kernel_MSD_GPU_final_regular(MSD_conf->final_gridSize,MSD_conf->final_blockSize,
				    &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE],
				    d_MSD, MSD_conf->nBlocks_total);
#ifdef MSD_DEBUG
  cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
  printf("Before outlier rejection: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
  printf("---------------------------<\n");
#endif
  for(int i=0; i<5; i++){
    call_kernel_MSD_BLN_calculate_partials_2d_and_minmax_with_outlier_rejection(MSD_conf->partials_gridSize,MSD_conf->partials_blockSize, d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD,  MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset, OR_sigma_multiplier);
    call_kernel_MSD_GPU_final_nonregular(MSD_conf->final_gridSize,MSD_conf->final_blockSize,
					 &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE],
					 d_MSD, MSD_conf->nBlocks_total);
#ifdef MSD_DEBUG
    cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
    printf("Rejection %d: Mean: %e, Standard deviation: %e; Elements:%zu;\n", i, h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
    printf("---------------------------<\n");
#endif
  }
  call_kernel_MSD_GPU_final_nonregular(MSD_conf->final_gridSize,MSD_conf->final_blockSize,
				       &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE],
				       d_MSD, d_previous_partials, MSD_conf->nBlocks_total);
	
#ifdef MSD_DEBUG
  cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
  printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
  printf("---------------------------<\n");
#endif
	
  return(0);
}


int MSD_outlier_rejection_continuous(float *d_MSD, float *d_input, float *d_previous_partials, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier) {
  int result;
  MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
  float *d_temp;
  cudaMalloc((void **) &d_temp, conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float));
  result = MSD_outlier_rejection_continuous(d_MSD, d_input, d_previous_partials, d_temp, &conf, OR_sigma_multiplier);
  cudaFree(d_temp);
  return(result);
}


//MSD_BLN_pw_continuous_OR
int MSD_outlier_rejection_grid(float *d_MSD, float *d_input, float *d_previous_partials, float *d_temp, MSD_Configuration *MSD_conf, float OR_sigma_multiplier){
#ifdef MSD_DEBUG
  float h_MSD[3];
  MSD_conf->print();
#endif
	
  MSD_init();
  call_kernel_MSD_GPU_calculate_partials_2d_and_minmax(MSD_conf->partials_gridSize,MSD_conf->partials_blockSize, d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int)MSD_conf->offset);
  call_kernel_MSD_GPU_final_regular(MSD_conf->final_gridSize,MSD_conf->final_blockSize,
				    &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE],
				    d_MSD, MSD_conf->nBlocks_total);
#ifdef MSD_DEBUG
  cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
  printf("Before outlier rejection: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
  printf("---------------------------<\n");
#endif
  for(int i=0; i<5; i++){
    call_kernel_MSD_BLN_calculate_partials_2d_and_minmax_with_outlier_rejection(MSD_conf->partials_gridSize,MSD_conf->partials_blockSize, d_input, &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE], d_MSD,  MSD_conf->nSteps.y, (int) MSD_conf->nTimesamples, (int) MSD_conf->offset, OR_sigma_multiplier);
    call_kernel_MSD_GPU_final_nonregular(MSD_conf->final_gridSize,MSD_conf->final_blockSize,
					 &d_temp[MSD_conf->address*MSD_PARTIAL_SIZE],
					 d_MSD, MSD_conf->nBlocks_total);
#ifdef MSD_DEBUG
    cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
    printf("Rejection %d: Mean: %e, Standard deviation: %e; Elements:%zu;\n", i, h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
    printf("---------------------------<\n");
#endif
  }
	
#ifdef MSD_DEBUG
  cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
  printf("Before grid rejection: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
  printf("---------------------------<\n");
#endif
	
  call_kernel_MSD_BLN_grid_outlier_rejection_new(MSD_conf->final_gridSize, MSD_conf->final_blockSize, d_temp, d_MSD, MSD_conf->nBlocks_total+MSD_conf->address, OR_sigma_multiplier);
	
#ifdef MSD_DEBUG
  cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
  printf("Output: Mean: %e, Standard deviation: %e; Elements:%zu;\n", h_MSD[0], h_MSD[1], (size_t) h_MSD[2]);
  printf("---------------------------<\n");
#endif
	
  return(0);
}

//------------- MSD with outlier rejection
//---------------------------------------------------------------




void Find_MSD(float *d_MSD, float *d_input, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, int enable_outlier_rejection){
  MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
  float *d_MSD_workarea;
  cudaError_t e = cudaMalloc((void **) &d_MSD_workarea, conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float));
  
  if(e != cudaSuccess) {
    LOG(log_level::error, "Could not cudaMemcpy in aa_device_MSD.cu (" + std::string(cudaGetErrorString(e)) + ")");
  }
  
  if(enable_outlier_rejection){
    MSD_outlier_rejection(d_MSD, d_input, d_MSD_workarea, &conf, OR_sigma_multiplier);
  }
  else {
    MSD_normal(d_MSD, d_input, d_MSD_workarea, &conf);
  }

  e = cudaFree(d_MSD_workarea);

  if(e != cudaSuccess) {
    LOG(log_level::error, "Could not cudaMemcpy in aa_device_MSD.cu (" + std::string(cudaGetErrorString(e)) + ")");
  }
  
}

void Find_MSD(float *d_MSD, float *d_input, float *d_MSD_workarea, MSD_Configuration *conf, float OR_sigma_multiplier, int enable_outlier_rejection){
  if(enable_outlier_rejection){
    MSD_outlier_rejection(d_MSD, d_input, d_MSD_workarea, conf, OR_sigma_multiplier);
  }
  else {
    MSD_normal(d_MSD, d_input, d_MSD_workarea, conf);
  }
}

void Find_MSD_continuous(float *d_MSD, float *d_input, float *d_previous_partials, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, int enable_outlier_rejection){
  MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
  float *d_MSD_workarea;
  cudaError_t e = cudaMalloc((void **) &d_MSD_workarea, conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float));

  if(e != cudaSuccess) {
    LOG(log_level::error, "Could not cudaMemcpy in aa_device_MSD.cu (" + std::string(cudaGetErrorString(e)) + ")");
  }
  
  if(enable_outlier_rejection){
    MSD_outlier_rejection_continuous(d_MSD, d_input, d_previous_partials, d_MSD_workarea, &conf, OR_sigma_multiplier);
  }
  else {
    MSD_normal_continuous(d_MSD, d_input, d_previous_partials, d_MSD_workarea, &conf);
  }

  e = cudaFree(d_MSD_workarea);
  
  if(e != cudaSuccess) {
    LOG(log_level::error, "Could not cudaMemcpy in aa_device_MSD.cu (" + std::string(cudaGetErrorString(e)) + ")");
  }
  
}

void Find_MSD_continuous(float *d_MSD, float *d_input, float *d_previous_partials, float *d_MSD_workarea, MSD_Configuration *conf, float OR_sigma_multiplier, int enable_outlier_rejection){
  if(enable_outlier_rejection){
    MSD_outlier_rejection_continuous(d_MSD, d_input, d_previous_partials, d_MSD_workarea, conf, OR_sigma_multiplier);
  }
  else {
    MSD_normal_continuous(d_MSD, d_input, d_previous_partials, d_MSD_workarea, conf);
  }
}



//---------------------------------------------------------------
//------------- MSD with outlier rejection on grid

int MSD_grid_outlier_rejection(float *d_MSD, float *d_input, int CellDim_x, int CellDim_y, int nTimesamples, int nDMs, int offset, float multiplier){
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
  call_kernel_MSD_BLN_grid_calculate_partials(gridSize,blockSize,nThreads*8, d_input, d_output, x_steps, y_steps, nTimesamples, 0);
  call_kernel_MSD_BLN_grid_outlier_rejection(final_gridSize, final_blockSize, d_output, d_MSD, GridSize_x*GridSize_y, (float) (CellDim_x*CellDim_y), multiplier);

  //---------> De-allocation of temporary memory
  cudaFree(d_output);
	
  return(1);
}

} //namespace astroaccelerate
