#include <iostream>

#include "aa_fdas_host.hpp"
#include "aa_log.hpp"

namespace astroaccelerate {

  /** \brief Print fdas parameters. */
  void  fdas_print_params_h()
  {
    printf("\n\nParameters defined in params.h:\n\t-------------------\n");
    //  printf("\nSampling time: TSAMP %g\n", TSAMP);
    printf("\nSpeed of light: SLIGHT %g\n", SLIGHT);
    printf("\nTemplate length for FFT: KERNLEN = RADIX*POTWO %d\n", KERNLEN);
    printf("\nAcceleration step in fourier bins (z): ACCEL_STEP %f\n", ACCEL_STEP);
    printf("\nAcceleration step in fourier bins (z) reciprocal: ACCEL_STEP_R %f\n", ACCEL_STEP_R);
    printf("\nMaximum acceleration in fourier bins (z): ZMAX %d\n", ZMAX);
    printf("\nNumber of templates including zero acceleration: NKERN %d\n", NKERN);
    //  printf("\nLowest acceleration in fourier bins (z) (for harmonic sum): ZLO %d\n", ZLO);
    printf("\nThread block size in x direction for 2-D thread block convolution GPU kernels : TBSIZEX %d\n", TBSIZEX);
    printf("\nThread block size in Y direction for 2-D thread block convolution GPU kernels : TBSIZEY %d\n", TBSIZEY);
    printf("\nThread block size in x direction for 2-D thread block power spectrum GPU kernels : PTBSIZEX %d\n", PTBSIZEX);
    printf("\nThread block size in y direction for 2-D thread block power spectrum GPU kernels : PTBSIZEY %d\n", PTBSIZEY);
    printf("\n\nCustom FFT specific parameters:\n\t------------------\n" );
    printf("\n\n\t--------------\n\n");
  }

  /** \brief Check CUDA devices. */
  void fdas_cuda_check_devices(int devid) {
    int devcount;
    cudaError_t e = cudaGetDeviceCount(&devcount);
    
    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaGetDeviceCount in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    printf("\nDetected %d CUDA Capable device(s)\n", devcount);
  }

  /** \brief Allocate GPU arrays for fdas. */
  void fdas_alloc_gpu_arrays(fdas_gpuarrays *arrays,  cmd_args *cmdargs)
  {
    printf("\nAllocating gpu arrays:\n"); 

    if (cmdargs->inbin){
      printf("\nF-fdot array will be interbinned\n");
    }
    double gbyte = 1024.0*1024.0*1024.0;
    //double mbyte = 1024.0*1024.0;

    // Memory allocations for gpu real fft input / output signal
    cudaError_t e = cudaMalloc((void**)&arrays->d_in_signal, arrays->mem_insig);
    
    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMalloc in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }

    e = cudaMalloc((void**)&arrays->d_fft_signal, arrays->mem_rfft);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMalloc in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }

    //Allocating arrays for fourier domain convolution  
    e = cudaMalloc((void**)&arrays->d_ext_data, arrays->mem_extsig);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMalloc in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }

    //templates
    e = cudaMalloc((void**)&arrays->d_kernel, KERNLEN*sizeof(float2)*NKERN );

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMalloc in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }

    //ffdot planes
    e = cudaMalloc((void**)&arrays->d_ffdot_pwr, arrays->mem_ffdot );

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMalloc in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    //initialise array
    e = cudaMemset(arrays->d_ffdot_pwr, 0, arrays->mem_ffdot);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemset in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }

    printf("ffdot x size: %lu",(unsigned long)arrays->mem_ffdot/sizeof(float)/(unsigned long)NKERN);
    if(cmdargs->basic==1){
      e = cudaMalloc(&arrays->d_ffdot_cpx, arrays->mem_ffdot_cpx);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
    }

    if(cmdargs->kfft && cmdargs->inbin){
      //    printf("mem_ipedge = %u ",mem_ipedge/);
      e = cudaMalloc(&arrays->ip_edge_points, arrays->mem_ipedge);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMalloc in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
    }
   
    // Added by KA
    if ( cudaSuccess != cudaMalloc((void**) &arrays->d_fdas_peak_list, arrays->mem_max_list_size)) printf("Allocation error in FDAS: d_fdas_peak_list\n");
	
    // check allocated/free memory
    size_t mfree,  mtotal;
    e = cudaMemGetInfo ( &mfree, &mtotal );

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemGetInfo in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    printf("\nMemory allocation finished: Total memory for this device: %.2f GB\nAvailable memory left on this device: %.2f GB \n", mtotal/gbyte, mfree/gbyte);
  }

  /** \brief Free GPU arrays for fdas. */
  void fdas_free_gpu_arrays(fdas_gpuarrays *arrays,  cmd_args *cmdargs)
  {

    cudaError_t e = cudaFree(arrays->d_in_signal);
    
    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaFree in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    e = cudaFree(arrays->d_fft_signal);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaFree in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    e = cudaFree(arrays->d_ext_data);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaFree in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    e = cudaFree(arrays->d_ffdot_pwr);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaFree in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    e = cudaFree(arrays->d_kernel);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaFree in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    if(cmdargs->basic) {
      e = cudaFree(arrays->d_ffdot_cpx);
      
      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaFree in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
      } 
    }

    if(cmdargs->kfft && cmdargs->inbin) {
      e = cudaFree(arrays->ip_edge_points);
      
      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaFree in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
      } 
    }
	
    // Added by KA
    cudaFree(arrays->d_fdas_peak_list);
  }

  /**
   * \brief Create kernel templates for the correlation technique (Ransom et. al. 2002), and upload + FFT to GPU memory.
   * \brief Using functions from the original PRESTO accelsearch code, (small adaptations for variables and remove normal interpolation management - input is already interpolated signal).
   * \author Scott Ransom.
   */
  void fdas_create_acc_kernels(cufftComplex* d_kernel, cmd_args *cmdargs ) {
    int ii;
    int inbin = 1;
    cufftComplex *h_kernel, *tempkern;
    cufftHandle templates_plan; // for host kernel fft
    int nrank = 1;
    int n[] = {KERNLEN};
    int idist = n[0], odist =n[0];
    int *inembed = n, *onembed = n;
    int istride =1, ostride = 1;

    //allocate kernel array and prepare fft
    h_kernel = (cufftComplex*) malloc(NKERN*KERNLEN*sizeof(float2));

    // batched fft plan for the templates array
    cufftPlanMany( &templates_plan, nrank, n, inembed , istride, 
		   idist, onembed, ostride,
		   odist, CUFFT_C2C, NKERN); 

    for (ii = 0; ii < NKERN; ii++){
      double z = (-ZMAX+ii*ACCEL_STEP);
      int halfwidth = presto_z_resp_halfwidth(z, LOWACC) ;
      int numkern = 2 * halfwidth * inbin;
      tempkern = presto_gen_z_response(0.0, inbin, z, numkern);
      presto_place_complex_kernel(tempkern, numkern, (h_kernel+ii*KERNLEN), KERNLEN);
      free(tempkern);
    }
  
    //!TEST!: replace templates here. Template width: numkern; padded width: KERNLEN
#ifdef FDAS_CONV_TEST
    for (ii = 0; ii < NKERN; ii++){
      int boxcar_width=ii*FDAS_TEST_FILTER_INCREMENT;
      for(int f=0; f<KERNLEN; f++){
	h_kernel[ii*KERNLEN + f].x = 0;
	h_kernel[ii*KERNLEN + f].y = 0;
	  
	if(f<boxcar_width/2) h_kernel[ii*KERNLEN + f].x = 1.0;
	if(f>=(KERNLEN-boxcar_width/2)) h_kernel[ii*KERNLEN + f].x = 1.0;
      }
    }
#endif
    //!TEST!: replace templates here. Template width: numkern; padded width: KERNLEN
  
    cudaError_t e = cudaMemcpy( d_kernel, h_kernel, KERNLEN*sizeof(float2)* NKERN, cudaMemcpyHostToDevice); // upload kernels to GPU
    
    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
#ifndef NOCUST
    //use kerel's non-reordered fft
    if (cmdargs->kfft)
      call_kernel_customfft_fwd_temps_no_reorder(d_kernel);  
#endif
    //use cuFFT to transform the templates
    if (cmdargs->basic)
      cufftExecC2C(templates_plan, d_kernel, d_kernel, CUFFT_FORWARD); 

    free(h_kernel);

  }

  /** \brief Create CUDA cufft fftplans for FDAS. */
  void fdas_cuda_create_fftplans(fdas_cufftplan *fftplans, fdas_params *params) {
    /*check plan memory overhead and create plans */
    double mbyte = 1024.0*1024.0;
    //double gbyte = mbyte*1024.0;
 
    //set cufft plan parameters
    size_t sig_worksize, real_worksize;
    int nrank = 1;
    int n[] = {KERNLEN};
    int idist = n[0], odist =n[0];
    int *inembed = n, *onembed = n;
    int istride =1, ostride = 1;

    //estimate plan memory for real fft
    cufftResult e = cufftEstimate1d( params->nsamps, CUFFT_R2C, 1, &real_worksize);
    
    if(e != CUFFT_SUCCESS) {
      LOG(log_level::error, "Could not cufftEstimate1d in aa_fdas_host.cu");
    }
    
    printf("\nsignal real fft plan requires extra %f MB of memory\n", real_worksize / mbyte);

    //estimate plan memory for forward fft
    e = cufftEstimateMany(nrank, n,inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, params->nblocks, &sig_worksize);

    if(e != CUFFT_SUCCESS) {
      LOG(log_level::error, "Could not cudaEstimateMany in aa_fdas_host.cu");
    }
    
    printf("\nsignal forward fft plan requires extra  %f MB of memory\n the same plan is used for the inverse fft", sig_worksize / mbyte);
  
    // real plan
    size_t rworksize;
    int rn[] = {params->nsamps};
    int *rinembed = rn, *ronembed = rn;
    int ridist = rn[0], rodist = params->rfftlen;
 
    cufftCreate(&fftplans->realplan);
    e = cufftMakePlanMany( fftplans->realplan, nrank, rn, rinembed, istride, ridist, ronembed, ostride, rodist, CUFFT_R2C, 1, &rworksize);

    if(e != CUFFT_SUCCESS) {
      LOG(log_level::error, "Could not cufftMakePlanMany in aa_fdas_host.cu");
    }
    
    cudaDeviceSynchronize();
    //getLastCudaError("\nCuda Error real fft plan\n");

    // forward batched plan - same used for inverse
    e = cufftCreate(&fftplans->forwardplan);

    if(e != CUFFT_SUCCESS) {
      LOG(log_level::error, "Could not cufftCreate in aa_fdas_host.cu");
    }
    
    e = cufftMakePlanMany( fftplans->forwardplan, nrank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, params->nblocks, &sig_worksize);

    if(e != CUFFT_SUCCESS) {
      LOG(log_level::error, "Could not cufftMakePlanMany in aa_fdas_host.cu");
    }
    
    cudaDeviceSynchronize();
    //getLastCudaError("\nCuda Error forward fft plan\n");
    printf("\ncuFFT plans done \n");
  }


  /** \brief Perform basic fourier domain accelerated search (fdas). */
  void fdas_cuda_basic(fdas_cufftplan *fftplans, fdas_gpuarrays *gpuarrays, cmd_args *cmdargs, fdas_params *params)
  {
    /* Basic GPU fdas algorithm using cuFFT */
    //int inbin;
    int cthreads = TBSIZEX;
    int cblocks = KERNLEN/TBSIZEX;

    dim3 pwthreads(PTBSIZEX, PTBSIZEY);
    dim3 pwblocks((params->sigblock / PTBSIZEX) + 1, NKERN/PTBSIZEY);

    /* if (cmdargs->inbin)
       inbin = 2;
       else
       inbin = 1;
    */
    //real fft
#ifndef FDAS_CONV_TEST
    cufftExecR2C(fftplans->realplan, gpuarrays->d_in_signal, gpuarrays->d_fft_signal);
#endif
  
#ifdef FDAS_CONV_TEST
    float2 *f2temp;
    float *ftemp;
    ftemp  = (float *)malloc(params->rfftlen*sizeof(float));
    f2temp = (float2 *)malloc(params->rfftlen*sizeof(float2));
    cudaError_t e = cudaMemcpy(ftemp, gpuarrays->d_in_signal, (params->rfftlen)*sizeof(float), cudaMemcpyDeviceToHost);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    for(int f=0; f<params->rfftlen; f++){
      f2temp[f].x = ftemp[f];
      f2temp[f].y = 0;
    }
    e = cudaMemcpy(gpuarrays->d_fft_signal, f2temp, (params->rfftlen)*sizeof(float2), cudaMemcpyHostToDevice);
    
    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    free(ftemp);
    free(f2temp);
#endif

    if (cmdargs->norm){
      //  PRESTO deredden - remove red noise.
      // TODO: replace with GPU version
      float2 *fftsig;
      fftsig = (float2*)malloc((params->rfftlen)*sizeof(float2)); 
    
      cudaError_t e = cudaMemcpy(fftsig, gpuarrays->d_fft_signal, (params->rfftlen)*sizeof(float2), cudaMemcpyDeviceToHost);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
      
      presto_dered_sig(fftsig, params->rfftlen);
      e = cudaMemcpy(gpuarrays->d_fft_signal, fftsig, (params->rfftlen)*sizeof(float2), cudaMemcpyHostToDevice);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
      
      free(fftsig);
    }

    //overlap-copy
    call_kernel_cuda_overlap_copy(gpuarrays->d_ext_data, gpuarrays->d_fft_signal, params->sigblock, params->rfftlen, params->extlen, params->offset, params->nblocks );

    if (cmdargs->norm){
      //  PRESTO block median normalization
      // TODO: replace with GPU version
      float2 *extsig;
      extsig = (float2*)malloc((params->extlen)*sizeof(float2));
      cudaError_t e = cudaMemcpy(extsig, gpuarrays->d_ext_data, (params->extlen)*sizeof(float2), cudaMemcpyDeviceToHost);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
      
      for(int b=0; b<params->nblocks; ++b)
	presto_norm(extsig+b*KERNLEN, KERNLEN);
      e = cudaMemcpy(gpuarrays->d_ext_data, extsig, (params->extlen)*sizeof(float2), cudaMemcpyHostToDevice);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
      
      free(extsig);
    }

    //complex block fft
    cufftExecC2C(fftplans->forwardplan, gpuarrays->d_ext_data, gpuarrays->d_ext_data, CUFFT_FORWARD);

    //complex multiplication kernel
    call_kernel_cuda_convolve_reg_1d_halftemps(cblocks, cthreads, gpuarrays->d_kernel, gpuarrays->d_ext_data, gpuarrays->d_ffdot_cpx, params->extlen, params->scale);

    //inverse fft
    for (int k=0; k < ZMAX/2; k++){
      cufftExecC2C(fftplans->forwardplan, gpuarrays->d_ffdot_cpx + k * params->extlen, gpuarrays->d_ffdot_cpx + k *params->extlen, CUFFT_INVERSE);
      cufftExecC2C(fftplans->forwardplan, gpuarrays->d_ffdot_cpx + (ZMAX-k) * params->extlen, gpuarrays->d_ffdot_cpx + (ZMAX-k) *params->extlen, CUFFT_INVERSE);
    }
    // z=0
    cufftExecC2C(fftplans->forwardplan, gpuarrays->d_ffdot_cpx + ((ZMAX/2) * params->extlen), gpuarrays->d_ffdot_cpx + ((ZMAX/2) * params->extlen), CUFFT_INVERSE);

    //power spectrum 
    if (cmdargs->inbin){
      call_kernel_cuda_ffdotpow_concat_2d_inbin(pwblocks, pwthreads, gpuarrays->d_ffdot_cpx, gpuarrays->d_ffdot_pwr, params->sigblock, params->offset, params->nblocks, params->extlen, params->siglen);
    }
    else{
      call_kernel_cuda_ffdotpow_concat_2d(pwblocks, pwthreads, gpuarrays->d_ffdot_cpx, gpuarrays->d_ffdot_pwr, params->sigblock, params->offset, params->nblocks, params->extlen, params->siglen);
    }
  }

#ifndef NOCUST
  void fdas_cuda_customfft(fdas_cufftplan *fftplans, fdas_gpuarrays *gpuarrays, cmd_args *cmdargs, fdas_params *params) {
    //int nthreads;
    dim3 cblocks(params->nblocks, NKERN/2); 

    //real fft
#ifndef FDAS_CONV_TEST
    cufftExecR2C(fftplans->realplan, gpuarrays->d_in_signal, gpuarrays->d_fft_signal);
#endif

#ifdef FDAS_CONV_TEST
    float2 *f2temp;
    float *ftemp;
    ftemp  = (float *)malloc(params->rfftlen*sizeof(float));
    f2temp = (float2 *)malloc(params->rfftlen*sizeof(float2));
    cudaError_t e = cudaMemcpy(ftemp, gpuarrays->d_in_signal, (params->rfftlen)*sizeof(float), cudaMemcpyDeviceToHost);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    for(int f=0; f<params->rfftlen; f++){
      f2temp[f].x = ftemp[f];
      f2temp[f].y = 0;
    }
    e = cudaMemcpy(gpuarrays->d_fft_signal, f2temp, (params->rfftlen)*sizeof(float2), cudaMemcpyHostToDevice);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    free(ftemp);
    free(f2temp);
#endif
  

    if (cmdargs->norm){
      //  PRESTO deredden - remove red noise.
      // TODO: replace with GPU version
      float2 *fftsig;
      fftsig = (float2*)malloc((params->rfftlen)*sizeof(float2)); 
    
      cudaError_t e = cudaMemcpy(fftsig, gpuarrays->d_fft_signal, (params->rfftlen)*sizeof(float2), cudaMemcpyDeviceToHost);
      
      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
      
      presto_dered_sig(fftsig, params->rfftlen);
      e = cudaMemcpy(gpuarrays->d_fft_signal, fftsig, (params->rfftlen)*sizeof(float2), cudaMemcpyHostToDevice);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
      
      free(fftsig);
    }

    //overlap-copy
    call_kernel_cuda_overlap_copy_smallblk(params->nblocks, gpuarrays->d_ext_data, gpuarrays->d_fft_signal, params->sigblock, params->rfftlen, params->extlen, params->offset, params->nblocks );

    if (cmdargs->norm){
      //  PRESTO block median normalization
      // TODO: replace with GPU version
      float2 *extsig;
      extsig = (float2*)malloc((params->extlen)*sizeof(float2));
      cudaError_t e = cudaMemcpy(extsig, gpuarrays->d_ext_data, (params->extlen)*sizeof(float2), cudaMemcpyDeviceToHost);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
      
      for(int b=0; b<params->nblocks; ++b)
	presto_norm(extsig+b*KERNLEN, KERNLEN);
      e = cudaMemcpy(gpuarrays->d_ext_data, extsig, (params->extlen)*sizeof(float2), cudaMemcpyHostToDevice);

      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
      
      free(extsig);
    }

    // Custom FFT convolution kernel
    if(cmdargs->inbin){
      call_kernel_cuda_convolve_customfft_wes_no_reorder02_inbin(params->nblocks, gpuarrays->d_kernel, gpuarrays->d_ext_data, gpuarrays->d_ffdot_pwr, params->sigblock, params->extlen, params->siglen, params->offset, params->scale, gpuarrays->ip_edge_points);
    }
    else{
      //cuda_convolve_customfft_wes_no_reorder02<<< params->nblocks, KERNLEN >>>( gpuarrays->d_kernel, gpuarrays->d_ext_data, gpuarrays->d_ffdot_pwr, params->sigblock, params->extlen, params->siglen, params->offset, params->scale);
		
      //-------------------------------------------
      dim3 gridSize(1, 1, 1);
      dim3 blockSize(1, 1, 1);
		
      /*
      //-------------------------------------------
      //Two elements per thread
      gridSize.x = params->nblocks;
      gridSize.y = 1;
      gridSize.z = 1;
      blockSize.x = KERNLEN/2;
      GPU_CONV_kFFT_mk11_2elem_2v<<<gridSize,blockSize>>>(gpuarrays->d_ext_data, gpuarrays->d_ffdot_pwr, gpuarrays->d_kernel, params->sigblock, params->offset, params->nblocks, params->scale);
      */
		
      //-------------------------------------------
      //Four elements per thread
      gridSize.x = params->nblocks;
      gridSize.y = 1;
      gridSize.z = 1;
      blockSize.x = KERNLEN/4;
      call_kernel_GPU_CONV_kFFT_mk11_4elem_2v(gridSize,blockSize, gpuarrays->d_ext_data, gpuarrays->d_ffdot_pwr, gpuarrays->d_kernel, params->sigblock, params->offset, params->nblocks, params->scale);
    }
  }
#endif

  /** \brief Write fdas list to disk. */
  void fdas_write_list(fdas_gpuarrays *gpuarrays, cmd_args *cmdargs, fdas_params *params, float *h_MSD, float dm_low, int dm_count, float dm_step, unsigned int list_size){
    int ibin=1;
    if (cmdargs->inbin) ibin=2;
    double tobs = (double)params->tsamp* (double)params->nsamps*ibin;
	
    if( !isnan(h_MSD[0]) || !isinf(h_MSD[0]) || !isnan(h_MSD[1]) || !isinf(h_MSD[1]) ){
      printf("Number of peaks:%d; mean:%f; strdev:%f\n", list_size, h_MSD[0], h_MSD[1]);
		
      float *h_fdas_peak_list = (float*)malloc(list_size*4*sizeof(float));
      cudaError_t e = cudaMemcpy(h_fdas_peak_list, gpuarrays->d_fdas_peak_list, list_size*4*sizeof(float), cudaMemcpyDeviceToHost);
      
      if(e != cudaSuccess) {
	LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
		
      //prepare file
      const char *dirname= "output_data";
      struct stat st = {0};

      if (stat(dirname, &st) == -1) {
	printf("\nDirectory %s does not exist, creating...\n", dirname);
	mkdir(dirname, 0700);
      }	
		
      FILE *fp_c;
      char pfname[200];
      sprintf(pfname, "acc_list_%f.dat", dm_low + ((float)dm_count)*dm_step);
      if ((fp_c=fopen(pfname, "w")) == NULL) {
	fprintf(stderr, "Error opening %s file for writing: %s\n",pfname, strerror(errno));
      }

      int i_list_size = (int)list_size;
      for(int f=0; f<i_list_size; f++){
	int j;
	double a, acc, acc1, jfreq, pow, SNR;
	a   = h_fdas_peak_list[4*f];
	j   = (int) h_fdas_peak_list[4*f + 1];
	pow = h_fdas_peak_list[4*f + 2];
	SNR = (pow-h_MSD[0])/h_MSD[1];
	jfreq = (double)(j) / tobs;
	acc = (double) (ZMAX - a* ACCEL_STEP);
	acc1 = acc*SLIGHT / jfreq / tobs / tobs;
	fprintf(fp_c, "%.2f\t%.3f\t%u\t%.3f\t%.3f\t%.3f\n", acc, acc1, j , jfreq, pow, SNR);
      }

      fclose(fp_c);
		
      free(h_fdas_peak_list);
    }
    else {
      printf("Error: mean or standard deviation was NaN or Inf!\n");
    }
  }


  /** \brief Write ffdot output data to disk. */
  void fdas_write_ffdot(fdas_gpuarrays *gpuarrays, cmd_args *cmdargs, fdas_params *params, float dm_low, int dm_count, float dm_step ) {
    int ibin=1;
    if (cmdargs->inbin)
      ibin=2;
    // Download, threshold and write ffdot data to file
    //int nsamps = params->nsamps;

    printf("\n\nWrite data for signal with %d samples\nf-fdot size=%u\n",params->nsamps, params->ffdotlen);
    float *h_ffdotpwr = (float*)malloc(params->ffdotlen* sizeof(float));
    //download data
    cudaError_t e = cudaMemcpy(h_ffdotpwr, gpuarrays->d_ffdot_pwr, params->ffdotlen*sizeof(float), cudaMemcpyDeviceToHost);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    // calculating statistics
    double total = 0.0;
    double mean;
    double stddev;
    // unsigned int j;
    int i_params_ffdotlen = (int)params->ffdotlen;
    for ( int j = 0; j < i_params_ffdotlen; ++j){
      total += (double)(h_ffdotpwr[j]);
      if(isnan(total)){
	printf("\nnan detected during sum for mean at j=%d\nValue at j:%f\n",j,h_ffdotpwr[j]);
	exit(1);
      }
    }
  
    mean = total / ((double)(params->ffdotlen)); 

    printf("\ntotal ffdot:%lf\tmean ffdot: %lf", total, mean);
      
    // Calculate standard deviation
    total = 0.0;
    for ( int j = 0; j < i_params_ffdotlen; ++j){
      total += ((double)h_ffdotpwr[j] - mean ) * ((double)h_ffdotpwr[j] - mean);
      if(isnan(total)||isinf(total)){
	printf("\ninf/nan detected during sum for mean at j=%d\nValue at j:%f\n",j,h_ffdotpwr[j]);
	exit(1);
      }
    }
    stddev = sqrt(abs(total) / (double)(params->ffdotlen - 1)); 
    printf("\nmean ffdot: %f\tstd ffdot: %lf\n", mean, stddev);

    //prepare file
    const char *dirname= "output_data";
    struct stat st = {0};

    if (stat(dirname, &st) == -1) {
      printf("\nDirectory %s does not exist, creating...\n", dirname);
      mkdir(dirname, 0700);
    }

    FILE *fp_c;
    char pfname[200];
    //  char *infilename;
    //  infilename = basename(cmdargs->afname);
    // filename needs to be acc_dm_%f, dm_low[i] + ((float)dm_count)*dm_step[i]
    //sprintf(pfname, "%s/out_inbin%d_%s",dirname,ibin,infilename);
    sprintf(pfname, "acc_%f.dat", dm_low + ((float)dm_count)*dm_step);
    printf("\nwriting results to file %s\n",pfname);
    if ((fp_c=fopen(pfname, "w")) == NULL) {
      fprintf(stderr, "Error opening %s file for writing: %s\n",pfname, strerror(errno));
      exit(1);
    }
    float pow, sigma;
    double tobs = (double)params->tsamp * (double)params->nsamps*ibin;
    unsigned int numindep = params->siglen*(NKERN+1)*ACCEL_STEP/6.95; // taken from PRESTO

    //write to file
    printf("\nWriting ffdot data to file...\n");

    for(int a = 0; a < NKERN; a++) {
      double acc = (double) (ZMAX - a* ACCEL_STEP);
      for( int j = 0; j < ibin*params->siglen; j++){
	pow =  h_ffdotpwr[a * ibin*params->siglen + j]; //(h_ffdotpwr[a * params->siglen + j]-mean)/stddev;
		
	if( pow > cmdargs->thresh) {
	  sigma = candidate_sigma(pow, cmdargs->nharms, numindep);//power, number of harmonics, number of independed searches=1...2^harms
	  //  sigma=1.0;
	  double jfreq = (double)(j) / tobs;
	  double acc1 = acc*SLIGHT / jfreq / tobs / tobs;
	  fprintf(fp_c, "%.2f\t%.3f\t%u\t%.3f\t%.3f\t%.3f\n", acc, acc1, j , jfreq, pow, sigma);
	}    
      }
    }

    fclose(fp_c);
    printf("\nFinished writing file %s\n",pfname);
    
    free(h_ffdotpwr);

  }


  /** \brief Write test ffdot to disk. */
  void fdas_write_test_ffdot(fdas_gpuarrays *gpuarrays, cmd_args *cmdargs, fdas_params *params, float dm_low, int dm_count, float dm_step ) {
    int ibin=1;
    if (cmdargs->inbin)
      ibin=2;
    /* Download, threshold and write ffdot data to file */
    //int nsamps = params->nsamps;

    printf("\n\nWrite data for signal with %d samples\nf-fdot size=%u\n",params->nsamps, params->ffdotlen);
    float *h_ffdotpwr = (float*)malloc(params->ffdotlen* sizeof(float));
    //download data
    cudaError_t e = cudaMemcpy(h_ffdotpwr, gpuarrays->d_ffdot_pwr, params->ffdotlen*sizeof(float), cudaMemcpyDeviceToHost);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpy in aa_fdas_host.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }

    // calculating statistics
    double total = 0.0;
    double mean;
    double stddev;
    // unsigned int j;
    int i_params_ffdotlen = (int)params->ffdotlen;
    for ( int j = 0; j < i_params_ffdotlen; ++j){
      total += (double)(h_ffdotpwr[j]);
      if(isnan(total)){
	printf("\nnan detected during sum for mean at j=%d\nValue at j:%f\n",j,h_ffdotpwr[j]);
	exit(1);
      }
    }
  
    mean = total / ((double)(i_params_ffdotlen)); 

    printf("\ntotal ffdot:%lf\tmean ffdot: %lf", total, mean);
      
    // Calculate standard deviation
    total = 0.0;
    for ( int j = 0; j < i_params_ffdotlen; ++j){
      total += ((double)h_ffdotpwr[j] - mean ) * ((double)h_ffdotpwr[j] - mean);
      if(isnan(total)||isinf(total)){
	printf("\ninf/nan detected during sum for mean at j=%d\nValue at j:%f\n",j,h_ffdotpwr[j]);
	exit(1);
      }
    }
    stddev = sqrt(abs(total) / (double)(i_params_ffdotlen - 1)); 
    printf("\nmean ffdot: %f\tstd ffdot: %lf\n", mean, stddev);

    //prepare file
    const char *dirname= "output_data";
    struct stat st = {0};

    if (stat(dirname, &st) == -1) {
      printf("\nDirectory %s does not exist, creating...\n", dirname);
      mkdir(dirname, 0700);
    }

    FILE *fp_c;
    char pfname[200];
    sprintf(pfname, "acc_fdas_conv_test.dat");
    printf("\nwriting results to file %s\n",pfname);
    if ((fp_c=fopen(pfname, "w")) == NULL) {
      fprintf(stderr, "Error opening %s file for writing: %s\n",pfname, strerror(errno));
      exit(1);
    }
    float pow;

    //write to file
    printf("\nWriting ffdot data to file...\n");

    for(int a = 0; a < NKERN; a++) {
      for( int j = 0; j < ibin*params->siglen; j++){
	pow =  h_ffdotpwr[a * ibin*params->siglen + j]; //(h_ffdotpwr[a * params->siglen + j]-mean)/stddev;
	fprintf(fp_c, "%u\t%u\t%f\n", a, j, pow); 
      }
    }

    fclose(fp_c);
    printf("\nFinished writing file %s\n",pfname);
    
    free(h_ffdotpwr);

  }
} //namespace astroaccelerate
