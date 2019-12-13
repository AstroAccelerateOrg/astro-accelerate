#include "aa_corner_turn.hpp"
#include "aa_gpu_timer.hpp"

namespace astroaccelerate {

void corner_turn(unsigned short *const d_input, float *const d_output, const int nchans, const int nsamp) {
    
    //{{{ Simple corner turn on the GPU
    
    int divisions_in_t = CT;
    int divisions_in_f = CF;
    int num_blocks_t = nsamp / divisions_in_t;
    int num_blocks_f = nchans / divisions_in_f;
    
    printf("\nCORNER TURN!");
    printf("\n%d %d", nsamp, nchans);
    printf("\n%d %d", divisions_in_t, divisions_in_f);
    printf("\n%d %d", num_blocks_t, num_blocks_f);
    
    dim3 threads_per_block(divisions_in_t, divisions_in_f);
    dim3 num_blocks(num_blocks_t, num_blocks_f);
    
    aa_gpu_timer timer;
   
	timer.Start();	 
    call_kernel_simple_corner_turn_kernel(num_blocks, threads_per_block, d_input, d_output, nchans, nsamp);
    cudaDeviceSynchronize();
    call_kernel_swap(num_blocks, threads_per_block, d_input, d_output, nchans, nsamp);
    cudaDeviceSynchronize();
	timer.Stop();
    
    double time = (double)timer.Elapsed();
    printf("\nPerformed CT: %lf (GPU estimate)", time/1000.0);
    printf("\nCT Gops based on %.2f ops per channel per tsamp: %lf", 10.0,(( 10.0*( divisions_in_t * divisions_in_f*num_blocks_t * num_blocks_f ) ) / ( time/1000.0 ) ) / 1000000000.0);
    printf("\nCT Device memory bandwidth in GB/s: %lf", ( ( sizeof(float) + sizeof(unsigned short) ) * ( divisions_in_t * divisions_in_f * num_blocks_t * num_blocks_f ) ) / ( time/1000.0 ) / 1000000000.0);
    
    //cudaMemcpy(d_input, d_output, inputsize, cudaMemcpyDeviceToDevice);
    //}}}
    
}


int corner_turn(float *const d_input, float *const d_output, const int primary_size, const int secondary_size) {
    int divisions_in_primary = 32;
    int divisions_in_secondary = 8;
    int num_blocks_primary = primary_size / divisions_in_primary;
    int num_blocks_secondary = secondary_size / divisions_in_secondary;
    
#ifdef CORNER_TURN_DEBUG
    printf("\nCORNER TURN!");
    printf("\n%d %d", primary_size, secondary_size);
    printf("\n%d %d", divisions_in_primary, divisions_in_secondary);
    printf("\n%d %d", num_blocks_primary, num_blocks_secondary);
#endif
    
    dim3 threads_per_block(divisions_in_primary, divisions_in_secondary);
    dim3 num_blocks(num_blocks_primary, num_blocks_secondary);
    call_kernel_simple_corner_turn_kernel(num_blocks, threads_per_block, d_input, d_output, primary_size, secondary_size);
    cudaDeviceSynchronize();
    
    return(0);
}

int corner_turn_SM(float *const d_input, float *const d_output, const int primary_size, const int secondary_size) {
    //---------> Task specific
    int nBlocks_x, nBlocks_y, nRest, Elements_per_block;
    
    Elements_per_block=CT_CORNER_BLOCKS*WARP;
    nBlocks_x=primary_size/Elements_per_block;
    nRest=primary_size - nBlocks_x*Elements_per_block;
    if(nRest>0) nBlocks_x++;
    
    nBlocks_y = secondary_size/WARP;
    if( secondary_size%WARP != 0 ) nBlocks_y++;
    
    //---------> CUDA block and CUDA grid parameters
    dim3 gridSize(nBlocks_x, nBlocks_y, 1);
    dim3 blockSize(CT_NTHREADS, 1, 1);
    
#ifdef CORNER_TURN_DEBUG
    int SM_size=(WARP*(WARP+1)*(CT_CORNER_BLOCKS+1))*4;
    printf("\n-------------- CORNER TURN DEBUG (corner_turn_SM) ----------------\n");
    printf("Grid  settings: x:%d; y:%d; z:%d;\n", gridSize.x, gridSize.y, gridSize.z);
    printf("Block settings: x:%d; y:%d; z:%d;\n", blockSize.x, blockSize.y, blockSize.z);
    printf("nRest:%d;\n", nRest);
    printf("Shared memory: %d bytes; %d floats\n", SM_size, SM_size/4);
#endif
    
    cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
    call_kernel_corner_turn_SM_kernel(gridSize, blockSize, d_input, d_output, primary_size, secondary_size);
    
    return(nRest);
}

} //namespace astroaccelerate
