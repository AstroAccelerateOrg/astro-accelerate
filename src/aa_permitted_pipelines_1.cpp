//
//  aa_permitted_pipelines_1.cpp
//  aapipeline
//
//  Created by Cees Carels on Friday 02/11/2018.
//  Copyright © 2018 Astro-Accelerate. All rights reserved.
//

#include "aa_permitted_pipelines_1.hpp"

void allocate_memory_gpu(const int &maxshift, const int &max_ndms, const int &nchans, int **const t_processed, unsigned short **const d_input, float **const d_output) {
    
    int time_samps = t_processed[0][0] + maxshift;
    printf("\n\n\n%d\n\n\n", time_samps);
    size_t gpu_inputsize = (size_t) time_samps * (size_t) nchans * sizeof(unsigned short);

    checkCudaErrors( cudaMalloc((void **) d_input, gpu_inputsize) );

    size_t gpu_outputsize = 0;
    if (nchans < max_ndms) {
        gpu_outputsize = (size_t)time_samps * (size_t)max_ndms * sizeof(float);
    }
    else {
        gpu_outputsize = (size_t)time_samps * (size_t)nchans * sizeof(float);
    }

    checkCudaErrors( cudaMalloc((void **) d_output, gpu_outputsize) );
    cudaMemset(*d_output, 0, gpu_outputsize);

}

void run_pipeline_1(const aa_filterbank_metadata &metadata, const aa_ddtr_strategy &ddtr_strategy, unsigned short *const input_buffer) {
    printf("NOTICE: Pipeline started run_pipeline_1.\n");
    int num_tchunks                     = ddtr_strategy.num_tchunks();
    int **t_processed                   = ddtr_strategy.t_processed();
    std::vector<float> dm_shifts        = ddtr_strategy.dmshifts();
    float* dmshifts                     = dm_shifts.data();
    int maxshift                        = ddtr_strategy.maxshift();
    int max_ndms                        = ddtr_strategy.max_ndms();
    int nchans                          = metadata.nchans();
    int nbits                           = metadata.nbits();
    int enable_zero_dm                  = 0;
    int enable_zero_dm_with_outliers    = 0;
    int failsafe                        = 0;
    long int inc                        = 0;
    float tsamp                         = metadata.tsamp();
    float tsamp_original                = tsamp;
    int maxshift_original               = maxshift;
    size_t range                        = ddtr_strategy.range();
    float tstart_local                  = 0.0;
    
    //Allocate GPU memory
    unsigned short *d_input             = NULL;
    float *d_output                     = NULL;

    allocate_memory_gpu(maxshift, max_ndms, nchans, t_processed, &d_input, &d_output);
    //Put the dm low, high, step struct contents into separate arrays again.
    //This is needed so that the kernel wrapper functions don't need to be modified.
    float *dm_low  = (float*)malloc(ddtr_strategy.range() * sizeof(float));
    float *dm_high = (float*)malloc(ddtr_strategy.range() * sizeof(float));
    float *dm_step = (float*)malloc(ddtr_strategy.range() * sizeof(float));
    int   *inBin   = (int*)malloc(ddtr_strategy.range() * sizeof(int));
    for(size_t i = 0; i < ddtr_strategy.range(); i++) {
        dm_low[i]   = ddtr_strategy.dm(i).low;
        dm_high[i]  = ddtr_strategy.dm(i).high;
        dm_step[i]  = ddtr_strategy.dm(i).step;
        inBin[i]    = ddtr_strategy.dm(i).inBin;
    }
    
    const int *ndms = ddtr_strategy.ndms_data();
    
    
    for(int t = 0; t < num_tchunks; t++) {
        printf("\nNOTICE: t_processed:\t%d, %d", t_processed[0][t], t);

        checkCudaErrors(cudaGetLastError());
        load_data(-1, inBin, d_input, &input_buffer[(long int) ( inc * nchans )], t_processed[0][t], maxshift, nchans, dmshifts);

        checkCudaErrors(cudaGetLastError());
        
        if (enable_zero_dm) {
            zero_dm(d_input, nchans, t_processed[0][t]+maxshift, nbits);
        }

        checkCudaErrors(cudaGetLastError());

        
        if (enable_zero_dm_with_outliers) {
            zero_dm_outliers(d_input, nchans, t_processed[0][t]+maxshift);
        }

        checkCudaErrors(cudaGetLastError());

        corner_turn(d_input, d_output, nchans, t_processed[0][t] + maxshift);
        
        checkCudaErrors(cudaGetLastError());

        int oldBin = 1;
        for(size_t dm_range = 0; dm_range < range; dm_range++) {
            printf("\n\nNOTICE: %f\t%f\t%f\t%d\n", ddtr_strategy.dm(dm_range).low, ddtr_strategy.dm(dm_range).high, ddtr_strategy.dm(dm_range).step, ddtr_strategy.ndms(dm_range));
            printf("\nAmount of telescope time processed: %f\n", tstart_local);
            maxshift = maxshift_original / inBin[dm_range];

            cudaDeviceSynchronize();
            checkCudaErrors(cudaGetLastError());

            load_data(dm_range, inBin, d_input, &input_buffer[(long int) ( inc * nchans )], t_processed[dm_range][t], maxshift, nchans, dmshifts);

            checkCudaErrors(cudaGetLastError());

            
            if (inBin[dm_range] > oldBin) {
                bin_gpu(d_input, d_output, nchans, t_processed[dm_range - 1][t] + maxshift * inBin[dm_range]);
                ( tsamp ) = ( tsamp ) * 2.0f;
            }
            
            checkCudaErrors(cudaGetLastError());

            dedisperse(dm_range, t_processed[dm_range][t], inBin, dmshifts, d_input, d_output, nchans, &tsamp, dm_low, dm_step, ndms, nbits, failsafe);

            checkCudaErrors(cudaGetLastError());

            oldBin = inBin[dm_range];
        }
        
        inc = inc + t_processed[0][t];
        printf("\nNOTICE: INC:\t%ld\n", inc);
        tstart_local = ( tsamp_original * inc );
        tsamp = tsamp_original;
        maxshift = maxshift_original;
    }
    printf("NOTICE: Pipeline ended run_pipeline_1.\n");
}
