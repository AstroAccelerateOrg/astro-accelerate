/* FDAS host functions */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include "headers/fdas_host.h"
#include "headers/params.h"
//#include <helper_functions.h>
#include <helper_cuda.h>
#include <curand.h>
#include <libgen.h>
//#include <random> // C++11 to use normal distribution

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
  printf("\nTAPS \t%d\n", TAPS);
  printf("\n\n\t--------------\n\n");
}

void fdas_cuda_check_devices(int devid)
{
  //int dev = 0;
  int devcount;
  //cudaDeviceProp deviceProp;

/* ******* Detect CUDA devices ******* */
  checkCudaErrors(cudaGetDeviceCount(&devcount));
  printf("\nDetected %d CUDA Capable device(s)\n", devcount);
/*
  for (dev = 0; dev < devcount; ++dev)
    {
      cudaSetDevice(dev);
      cudaGetDeviceProperties(&deviceProp, dev);
      printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);
    }
  if (devid<devcount){
    printf("\nsetting device %d (default)\n", devid);
    cudaSetDevice(devid);
  }
  else{
    printf("\nDevice %d not found, setting device 0 (default)\n", devid);
    cudaSetDevice(0);
  }
*/
//   cudaSetDevice(CARD);
}

void fdas_alloc_gpu_arrays(fdas_gpuarrays *arrays,  cmd_args *cmdargs)
{
  printf("\nAllocating gpu arrays:\n"); 

  if (cmdargs->inbin){
    printf("\nF-fdot array will be interbinned\n");
  }
    double gbyte = 1024.0*1024.0*1024.0;
    //double mbyte = 1024.0*1024.0;

  // Memory allocations for gpu real fft input / output signal
  checkCudaErrors(cudaMalloc((void**)&arrays->d_in_signal, arrays->mem_insig));

  checkCudaErrors(cudaMalloc((void**)&arrays->d_fft_signal, arrays->mem_rfft));

  //Allocating arrays for fourier domain convolution  
  checkCudaErrors(cudaMalloc((void**)&arrays->d_ext_data, arrays->mem_extsig)); 

  //templates
   checkCudaErrors(cudaMalloc((void**)&arrays->d_kernel, KERNLEN*sizeof(float2)*NKERN )); 

   //ffdot planes
   checkCudaErrors(cudaMalloc((void**)&arrays->d_ffdot_pwr, arrays->mem_ffdot ));
   //initialise array
   checkCudaErrors(cudaMemset(arrays->d_ffdot_pwr, 0, arrays->mem_ffdot));

   printf("ffdot x size: %zu",arrays->mem_ffdot/sizeof(float)/NKERN);
   if(cmdargs->basic==1){
     checkCudaErrors(cudaMalloc(&arrays->d_ffdot_cpx, arrays->mem_ffdot_cpx));
   }

   if(cmdargs->kfft && cmdargs->inbin){
         //    printf("mem_ipedge = %u ",mem_ipedge/);
     checkCudaErrors(cudaMalloc(&arrays->ip_edge_points, arrays->mem_ipedge));
   }
   
	// Added by KA
	if ( cudaSuccess != cudaMalloc((void**) &arrays->d_fdas_peak_list, arrays->mem_max_list_size)) printf("Allocation error in FDAS: d_fdas_peak_list\n");
	
   // check allocated/free memory
   size_t mfree,  mtotal;
   checkCudaErrors(cudaMemGetInfo ( &mfree, &mtotal ));
   printf("\nMemory allocation finished: Total memory for this device: %.2f GB\nAvailable memory left on this device: %.2f GB \n", mtotal/gbyte, mfree/gbyte);
}

void fdas_free_gpu_arrays(fdas_gpuarrays *arrays,  cmd_args *cmdargs)
{

    checkCudaErrors(cudaFree(arrays->d_in_signal));
    checkCudaErrors(cudaFree(arrays->d_fft_signal));
    checkCudaErrors(cudaFree(arrays->d_ext_data));
    checkCudaErrors(cudaFree(arrays->d_ffdot_pwr));
    checkCudaErrors(cudaFree(arrays->d_kernel));
    if(cmdargs->basic)
      checkCudaErrors(cudaFree(arrays->d_ffdot_cpx));

    if(cmdargs->kfft && cmdargs->inbin)
      checkCudaErrors(cudaFree(arrays->ip_edge_points));
	
	// Added by KA
	cudaFree(arrays->d_fdas_peak_list);
}
/*
void fdas_create_acc_sig(fdas_new_acc_sig *acc_sig, cmd_args *cmdargs)
/* Create accelerated signal with given parameters in a float array */
/*{
  double t0, tau;
  double omega = 2*M_PI*acc_sig->freq0;
  double accel;
  double tobs;
  // gaussian distribution from C++ <random>
  std::default_random_engine rgen;
  std::normal_distribution<float> gdist(0.0,cmdargs->nsig);
  tobs = (double) (TSAMP*acc_sig->nsamps);
  accel = ((double)acc_sig->zval * SLIGHT) / (acc_sig->freq0*tobs*tobs);
  printf("\n\npreparing test signal, observation time = %f s, %d nsamps f0 = %f Hz with %d harmonics\n", tobs, acc_sig->nsamps, acc_sig->freq0, acc_sig->nharms);
  printf("\nz = %d accelereation = %f m/s^2\n", acc_sig->zval, accel);
  acc_sig->acc_signal = (float*)malloc(acc_sig->nsamps*sizeof(float));
     
  printf("\nNow creating accelerated signal with fc=%f, accel=%f, harmonics=%d, duty cycle=%.1f%, noise=%d signal samples=%d, signal level: %.2f\n", acc_sig->freq0, accel, acc_sig->nharms, acc_sig->duty*100.0, cmdargs->nsig, acc_sig->nsamps,acc_sig->sigamp);
  
  for ( int i=0; i<acc_sig->nsamps; ++i){	    
    t0 = i*TSAMP;
    tau = t0 + (t0*(accel*t0) / SLIGHT /2.0);
    if (cmdargs->nsig!=0){
      acc_sig->acc_signal[i] = gdist(rgen);
    }
    for (int j = 1; j <= acc_sig->nharms; ++j){
      acc_sig->acc_signal[i] += (2.0/(j*M_PI)*sin(j*M_PI*acc_sig->duty))*acc_sig->sigamp*cos(j*omega*tau); 
    }
  }
  //Write to file 
  char afname[200];
  sprintf(afname, "data/acc_sig_8192x%d_%dharms_%dduty_%.3fHz_%dz_%dnsigma.dat",  acc_sig->mul,  acc_sig->nharms, (int)( acc_sig->duty*100.0),  acc_sig->freq0,  acc_sig->zval, acc_sig->nsig );
  write_output_file(afname, &acc_sig->acc_signal, acc_sig->nsamps );
  free(acc_sig->acc_signal);
}
*/

void fdas_create_acc_kernels(cufftComplex* d_kernel, cmd_args *cmdargs )
{
/* Create kernel templates for the correlation technique (Ransom et. al. 2002),  */
/* and upload + FFT to GPU memory. */
/* Using functions from the original PRESTO accelsearch code */
/* (small adaptations for variables and remove normal interpolation management  */
/* - input is already interpolated signal)   */
/* -credit to Scott Ransom  */
 
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
    tempkern = presto_gen_z_response( z, numkern, inbin);
    presto_place_complex_kernel(tempkern, numkern, (h_kernel+ii*KERNLEN), KERNLEN);
    free(tempkern);
  }
  checkCudaErrors( cudaMemcpy( d_kernel, h_kernel, KERNLEN*sizeof(float2)* NKERN, cudaMemcpyHostToDevice) ); // upload kernels to GPU

#ifndef NOCUST
  //use kerel's non-reordered fft
  if (cmdargs->kfft)
    customfft_fwd_temps_no_reorder<<<NKERN,KERNLEN>>>( d_kernel);  
#endif
  //use cuFFT to transform the templates
  if (cmdargs->basic)
    cufftExecC2C(templates_plan, d_kernel, d_kernel, CUFFT_FORWARD); 

  free(h_kernel);

}

void fdas_cuda_create_fftplans(fdas_cufftplan *fftplans, fdas_params *params)
{
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
  checkCudaErrors(cufftEstimate1d( params->nsamps, CUFFT_R2C, 1, &real_worksize));
  printf("\nsignal real fft plan requires extra %f MB of memory\n", real_worksize / mbyte);

  //estimate plan memory for forward fft
   checkCudaErrors(cufftEstimateMany(nrank, n,inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, params->nblocks, &sig_worksize));
  printf("\nsignal forward fft plan requires extra  %f MB of memory\n the same plan is used for the inverse fft", sig_worksize / mbyte);
  
  // real plan
  size_t rworksize;
  int rn[] = {params->nsamps};
  int *rinembed = rn, *ronembed = rn;
  int ridist = rn[0], rodist = params->rfftlen;
 
  cufftCreate(&fftplans->realplan);
  checkCudaErrors(cufftMakePlanMany( fftplans->realplan, nrank, rn, rinembed, istride, ridist, ronembed, ostride, rodist, CUFFT_R2C, 1, &rworksize));
  cudaDeviceSynchronize();
  getLastCudaError("\nCuda Error real fft plan\n");

  // forward batched plan - same used for inverse
  checkCudaErrors(cufftCreate(&fftplans->forwardplan));
  checkCudaErrors(cufftMakePlanMany( fftplans->forwardplan, nrank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, params->nblocks, &sig_worksize)); 
  cudaDeviceSynchronize();
  getLastCudaError("\nCuda Error forward fft plan\n");
  printf("\ncuFFT plans done \n");
}


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
  cufftExecR2C(fftplans->realplan, gpuarrays->d_in_signal, gpuarrays->d_fft_signal);

  if (cmdargs->norm){
    //  PRESTO deredden - remove red noise.
    // TODO: replace with GPU version
    float2 *fftsig;
    fftsig = (float2*)malloc((params->rfftlen)*sizeof(float2)); 
    
    checkCudaErrors( cudaMemcpy(fftsig, gpuarrays->d_fft_signal, (params->rfftlen)*sizeof(float2), cudaMemcpyDeviceToHost));
    presto_dered_sig(fftsig, params->rfftlen);
    checkCudaErrors( cudaMemcpy(gpuarrays->d_fft_signal, fftsig, (params->rfftlen)*sizeof(float2), cudaMemcpyHostToDevice));
    free(fftsig);
  }

  //overlap-copy
   cuda_overlap_copy<<<KERNLEN/64, 64 >>>(gpuarrays->d_ext_data, gpuarrays->d_fft_signal, params->sigblock, params->rfftlen, params->extlen, params->offset, params->nblocks );

   if (cmdargs->norm){
     //  PRESTO block median normalization
     // TODO: replace with GPU version
     float2 *extsig;
     extsig = (float2*)malloc((params->extlen)*sizeof(float2));
     checkCudaErrors( cudaMemcpy(extsig, gpuarrays->d_ext_data, (params->extlen)*sizeof(float2), cudaMemcpyDeviceToHost));
     for(int b=0; b<params->nblocks; ++b)
       presto_norm(extsig+b*KERNLEN, KERNLEN);
     checkCudaErrors( cudaMemcpy(gpuarrays->d_ext_data, extsig, (params->extlen)*sizeof(float2), cudaMemcpyHostToDevice));
     free(extsig);
   }

  //complex block fft
  cufftExecC2C(fftplans->forwardplan, gpuarrays->d_ext_data, gpuarrays->d_ext_data, CUFFT_FORWARD);

  //complex multiplication kernel

  cuda_convolve_reg_1d_halftemps<<<cblocks, cthreads >>>( gpuarrays->d_kernel, gpuarrays->d_ext_data, gpuarrays->d_ffdot_cpx, params->extlen, params->scale);

  //inverse fft
  for (int k=0; k < ZMAX/2; k++){
    cufftExecC2C(fftplans->forwardplan, gpuarrays->d_ffdot_cpx + k * params->extlen, gpuarrays->d_ffdot_cpx + k *params->extlen, CUFFT_INVERSE);
    cufftExecC2C(fftplans->forwardplan, gpuarrays->d_ffdot_cpx + (ZMAX-k) * params->extlen, gpuarrays->d_ffdot_cpx + (ZMAX-k) *params->extlen, CUFFT_INVERSE);
  }
  // z=0
  cufftExecC2C(fftplans->forwardplan, gpuarrays->d_ffdot_cpx + ((ZMAX/2) * params->extlen), gpuarrays->d_ffdot_cpx + ((ZMAX/2) * params->extlen), CUFFT_INVERSE);

  //power spectrum 
  if (cmdargs->inbin){
    cuda_ffdotpow_concat_2d_inbin<<< pwblocks, pwthreads >>>(gpuarrays->d_ffdot_cpx, gpuarrays->d_ffdot_pwr, params->sigblock, params->offset, params->nblocks, params->extlen, params->siglen);
  }
  else{
    cuda_ffdotpow_concat_2d <<< pwblocks, pwthreads >>>(gpuarrays->d_ffdot_cpx, gpuarrays->d_ffdot_pwr, params->sigblock, params->offset, params->nblocks, params->extlen, params->siglen);
  }
}

#ifndef NOCUST
void fdas_cuda_customfft(fdas_cufftplan *fftplans, fdas_gpuarrays *gpuarrays, cmd_args *cmdargs, fdas_params *params)
{
  //int nthreads;
  dim3 cblocks(params->nblocks, NKERN/2); 

  //real fft
  cufftExecR2C(fftplans->realplan, gpuarrays->d_in_signal, gpuarrays->d_fft_signal);

  if (cmdargs->norm){
    //  PRESTO deredden - remove red noise.
    // TODO: replace with GPU version
    float2 *fftsig;
    fftsig = (float2*)malloc((params->rfftlen)*sizeof(float2)); 
    
    checkCudaErrors( cudaMemcpy(fftsig, gpuarrays->d_fft_signal, (params->rfftlen)*sizeof(float2), cudaMemcpyDeviceToHost));
    presto_dered_sig(fftsig, params->rfftlen);
    checkCudaErrors( cudaMemcpy(gpuarrays->d_fft_signal, fftsig, (params->rfftlen)*sizeof(float2), cudaMemcpyHostToDevice));
    free(fftsig);
  }

  //overlap-copy
  cuda_overlap_copy_smallblk<<<params->nblocks, KERNLEN >>>(gpuarrays->d_ext_data, gpuarrays->d_fft_signal, params->sigblock, params->rfftlen, params->extlen, params->offset, params->nblocks );

  if (cmdargs->norm){
    //  PRESTO block median normalization
    // TODO: replace with GPU version
    float2 *extsig;
    extsig = (float2*)malloc((params->extlen)*sizeof(float2));
    checkCudaErrors( cudaMemcpy(extsig, gpuarrays->d_ext_data, (params->extlen)*sizeof(float2), cudaMemcpyDeviceToHost));
    for(int b=0; b<params->nblocks; ++b)
      presto_norm(extsig+b*KERNLEN, KERNLEN);
    checkCudaErrors( cudaMemcpy(gpuarrays->d_ext_data, extsig, (params->extlen)*sizeof(float2), cudaMemcpyHostToDevice));
    free(extsig);
  }

  // Custom FFT convolution kernel
  if(cmdargs->inbin){
    cuda_convolve_customfft_wes_no_reorder02_inbin<<< params->nblocks, KERNLEN >>>( gpuarrays->d_kernel, gpuarrays->d_ext_data, gpuarrays->d_ffdot_pwr, params->sigblock, params->extlen, params->siglen, params->offset, params->scale, gpuarrays->ip_edge_points);
  }
  else{
    cuda_convolve_customfft_wes_no_reorder02<<< params->nblocks, KERNLEN >>>( gpuarrays->d_kernel, gpuarrays->d_ext_data, gpuarrays->d_ffdot_pwr, params->sigblock, params->extlen, params->siglen, params->offset, params->scale);
  }
}
#endif

void fdas_write_list(fdas_gpuarrays *gpuarrays, cmd_args *cmdargs, fdas_params *params, float *h_MSD, float dm_low, int dm_count, float dm_step, unsigned int list_size, std::vector<float> &h_fdas_peak_list){
	int ibin=1;
	if (cmdargs->inbin) ibin=2;
	double tobs = (double)params->tsamp* (double)params->nsamps*ibin;
	
	if( !isnan(h_MSD[0]) || !isinf(h_MSD[0]) || !isnan(h_MSD[1]) || !isinf(h_MSD[1]) ){
		printf("Number of peaks:%d; mean:%f; strdev:%f\n", list_size, h_MSD[0], h_MSD[1]);
		
		//float *h_fdas_peak_list = (float*)malloc(list_size*4*sizeof(float));
		//checkCudaErrors(cudaMemcpy(h_fdas_peak_list, gpuarrays->d_fdas_peak_list, list_size*4*sizeof(float), cudaMemcpyDeviceToHost));

		unsigned current_size = h_fdas_peak_list.size();
		h_fdas_peak_list.resize(current_size+(list_size*4));


		float *h_fdas_peak_list = (float*)malloc(list_size*4*sizeof(float));
		checkCudaErrors(cudaMemcpy(&h_fdas_peak_list[current_size], gpuarrays->d_fdas_peak_list, list_size*4*sizeof(float), cudaMemcpyDeviceToHost));

		
		//prepare file
		const char *dirname= "output_data";
		struct stat st = {0};

		if (stat(dirname, &st) == -1) {
			printf("\nDirectory %s does not exist, creating...\n", dirname);
			mkdir(dirname, 0700);
		}	
		
		//FILE *fp_c;
		//char pfname[200];
		//sprintf(pfname, "acc_list_%f.dat", dm_low + ((float)dm_count)*dm_step);
		//if ((fp_c=fopen(pfname, "w")) == NULL) {
		//	fprintf(stderr, "Error opening %s file for writing: %s\n",pfname, strerror(errno));
		//}
		
		for(int f=0; f<list_size; f++){
			int j;
			double a, acc, acc1, jfreq, pow, SNR;
			a   = h_fdas_peak_list[current_size +(4*f)];
			j   = (int) h_fdas_peak_list[current_size + (4*f + 1)];
			pow = h_fdas_peak_list[current_size + (4*f + 2)];
			h_fdas_peak_list[current_size + (4*f + 1)] = (double)(j) / tobs;
			h_fdas_peak_list[current_size + (4*f + 2)] = (pow-h_MSD[0])/h_MSD[1];
			acc = (double) (ZMAX - a* ACCEL_STEP);
			acc1 = acc*SLIGHT / jfreq / tobs / tobs;
			h_fdas_peak_list[current_size + (4*f + 3)] = dm_low + ((float)dm_count)*dm_step;
			//fprintf(fp_c, "%.2f\t%.3f\t%u\t%.3f\t%.3f\t%.3f\n", acc, acc1, j , jfreq, pow, SNR);
		}

		//fclose(fp_c);
		
		//free(h_fdas_peak_list);
	}
	else {
		printf("Error: mean or standard deviation was NaN or Inf!\n");
	}
}

void fdas_write_ffdot(fdas_gpuarrays *gpuarrays, cmd_args *cmdargs, fdas_params *params, float dm_low, int dm_count, float dm_step ) {
  int ibin=1;
  if (cmdargs->inbin)
    ibin=2;
  /* Download, threshold and write ffdot data to file */
  //int nsamps = params->nsamps;

  printf("\n\nWrite data for signal with %d samples\nf-fdot size=%u\n",params->nsamps, params->ffdotlen);
  float *h_ffdotpwr = (float*)malloc(params->ffdotlen* sizeof(float));
  //download data
  checkCudaErrors(cudaMemcpy(h_ffdotpwr, gpuarrays->d_ffdot_pwr, params->ffdotlen*sizeof(float), cudaMemcpyDeviceToHost));

  // calculating statistics
  double total = 0.0;
  double mean;
  double stddev;
  // unsigned int j;
  for ( int j = 0; j < params->ffdotlen; ++j){
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
  for ( int j = 0; j < params->ffdotlen; ++j){
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
