//#define GPU_ANALYSIS_DEBUG

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "headers/params.h"

#include "headers/device_BC_plan.h"
#include "headers/device_peak_find.h"
#include "headers/device_MSD_Configuration.h"
#include "headers/device_MSD.h"
#include "headers/device_MSD_legacy.h"
#include "headers/device_SPS_long.h"
#include "headers/device_threshold.h"
#include "headers/device_single_FIR.h"

#include "timer.h"

struct MSD_Data {
	int taps;
	double mean;
	double sd;
};

//---------------------------------------------------------------------------------
//-------> Kahan MSD
void d_kahan_summation(float *signal, int nDMs, int nTimesamples, int offset, float *result, float *error){
	double sum;
	double sum_error;
	double a,b;
	
	sum=0;
	sum_error=0;
	for(int d=0;d<nDMs; d++){
		for(int s=0; s<(nTimesamples-offset); s++){
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

void d_kahan_sd(float *signal, int nDMs, int nTimesamples, int offset, double mean, float *result, float *error){
	double sum;
	double sum_error;
	double a,b,dtemp;
	
	sum=0;
	sum_error=0;
	for(int d=0;d<nDMs; d++){
		for(int s=0; s<(nTimesamples-offset); s++){
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


void MSD_Kahan(float *h_input, int nDMs, int nTimesamples, int offset, double *mean, double *sd){
	float error, signal_mean, signal_sd;
	int nElements=nDMs*(nTimesamples-offset);
	
	d_kahan_summation(h_input, nDMs, nTimesamples, offset, &signal_mean, &error);
	signal_mean=signal_mean/nElements;
	
	d_kahan_sd(h_input, nDMs, nTimesamples, offset, signal_mean, &signal_sd, &error);
	signal_sd=sqrt(signal_sd/nElements);

	*mean=signal_mean;
	*sd=signal_sd;
}

void MSD_on_GPU(float *h_input, float *d_input, float *d_MSD, float *signal_mean, float *signal_sd, float *signal_mean_bln, float *signal_sd_bln, float *signal_mean_bl_bln, float *signal_sd_bl_bln, int nDMs, int nTimesamples, int offset, float sigma_constant, float *MSD_limited_time, float *MSD_BLN_pw_time, float *MSD_BLN_grid_time){
	GpuTimer timer;
	float h_MSD[3];
	cudaMemcpy( d_input, h_input, ((size_t) nDMs*nTimesamples)*sizeof(float), cudaMemcpyHostToDevice);
	
	timer.Start();
	MSD_normal(d_input, d_MSD, nDMs, nTimesamples, offset);
	timer.Stop();
	(*MSD_limited_time) = timer.Elapsed();
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	(*signal_mean) = h_MSD[0];
	(*signal_sd)   = h_MSD[1];
	
	
	timer.Start();
	MSD_outlier_rejection(d_input, d_MSD, nDMs, nTimesamples, offset, sigma_constant);
	timer.Stop();
	(*MSD_BLN_pw_time) = timer.Elapsed();
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	(*signal_mean_bln) = h_MSD[0];
	(*signal_sd_bln)   = h_MSD[1];
	
	/*
	timer.Start();
	MSD_BLN_grid(d_input, d_MSD, 32, 32, nDMs, nTimesamples, offset, sigma_constant);
	timer.Stop();
	(*MSD_BLN_grid_time) = timer.Elapsed();
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	(*signal_mean_bl_bln) = h_MSD[0];
	(*signal_sd_bl_bln)   = h_MSD[1];
	*/
}

void MSD_on_GPU_LA(float *h_input, float *d_input, float *d_MSD, float *h_MSD_LA, float *h_MSD_BLN_LA, int nDMs, int nTimesamples, int offset, float sigma_constant){
	cudaMemcpy( d_input, h_input, ((size_t) nDMs*nTimesamples)*sizeof(float), cudaMemcpyHostToDevice);
	
	MSD_linear_approximation(d_input, d_MSD, 32, nDMs, nTimesamples, offset);
	cudaMemcpy(h_MSD_LA, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost);
	
	
	MSD_BLN_LA_pw_normal(d_input, d_MSD, 32, nDMs, nTimesamples, offset, sigma_constant);
	cudaMemcpy(h_MSD_BLN_LA, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
}

void MSD_on_GPU_halfed(float *h_input, float *d_input, float *d_MSD, float *signal_mean, float *signal_sd, float *signal_mean_bln, float *signal_sd_bln, int nDMs, int nTimesamples, int offset, float sigma_constant){
	float h_MSD[3];
	float *h_temp;
	int dt=nTimesamples/2;
	h_temp = new float[nDMs*dt];
	
	for(int d=0; d<nDMs; d++){
		for(int s=0; s<dt; s++){
			h_temp[d*dt + s] = h_input[d*nTimesamples + 2*s];
		}
	}
	
	cudaMemcpy( d_input, h_temp, ((size_t) nDMs*dt)*sizeof(float), cudaMemcpyHostToDevice);
	
	MSD_normal(d_input, d_MSD, nDMs, dt, offset/2);
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	(*signal_mean) = h_MSD[0];
	(*signal_sd)   = h_MSD[1];
	
	
	MSD_outlier_rejection(d_input, d_MSD, nDMs, dt, offset/2, sigma_constant);
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	(*signal_mean_bln) = h_MSD[0];
	(*signal_sd_bln)   = h_MSD[1];
	
	delete [] h_temp;
}
//-------> Kahan MSD
//---------------------------------------------------------------------------------


void Calculate_FIR(float *input, float *output, int nTaps, int nDMs, int nTimesamples, int ut) {
	int d,s,t;
	float ftemp;
	
	for(d=0; d<nDMs; d++){
		for(s=0; s<nTimesamples-ut; s++){
			ftemp=0;
			for(t=0; t<nTaps; t++){
				ftemp+=input[d*nTimesamples + s + t];
			}
			output[d*nTimesamples + s]=ftemp;
		}
	}	
}

void Decimate_in_time(float *h_input, float *h_CPU_decimate, int DIT_value, int DIT_factor, int nDMs, int nTimesamples, int offset){
	float ftemp;
	int decimated_timesamples;
	
	decimated_timesamples=nTimesamples/(DIT_value*DIT_factor);
	for(int d=0; d<nDMs; d++){
		for(int s=0; s<decimated_timesamples; s++){
			ftemp=0;
			for(int t=0; t<DIT_factor; t++){
				ftemp = ftemp + h_input[d*decimated_timesamples*DIT_factor + s*DIT_factor + t];
			}
			h_CPU_decimate[d*decimated_timesamples + s]=ftemp;
		}
	}
}

void Export_data(float *input, size_t nDMs, size_t nTimesamples, char *filename){
	FILE *fp_out;
	char mod_filename[200];
	
	sprintf(mod_filename,"%s.dat",filename);
	if (( fp_out = fopen(filename, "wb") ) == NULL)	{
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}
	fwrite(input, (nDMs*nTimesamples)*sizeof(float), 4, fp_out);
	fclose(fp_out);
	
	for(int d=0; d<nDMs; d++){
		sprintf(mod_filename,"%s_dm%d.dat",filename,d);
		if (( fp_out = fopen(filename, "wb") ) == NULL)	{
			fprintf(stderr, "Error opening output file!\n");
			exit(0);
		}
		fwrite(&input[d*nTimesamples], nTimesamples*sizeof(float), 4, fp_out);
		fclose(fp_out);		
	}
}

void export_file_nDM_nTimesamples(float *data, int nDMs, int nTimesamples, char *filename) {
	FILE *file_out;
	char str[200];

	sprintf(str, "%s_DM.dat", filename);
	if (( file_out = fopen(str, "w") ) == NULL)	{
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}

	printf("export nDMs\n");
	for (int s = 0; s < nTimesamples; s++) {
		for (int d = 0; d < nDMs; d++) {
			fprintf(file_out, "%f ", data[d*nTimesamples + s]);
		}
		fprintf(file_out, "\n");
	}

	fclose(file_out);

	sprintf(str, "%s_Time.dat", filename);
	if (( file_out = fopen(str, "w") ) == NULL)	{
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}

	printf("export nTimesamples\n");
	for (int d = 0; d < nDMs; d++) {
		for (int s = 0; s < nTimesamples; s++) {
			fprintf(file_out, "%f ", data[d*nTimesamples + s]);
		}
		fprintf(file_out, "\n");
	}

	fclose(file_out);
}


//---------------------------------------------------------------------------------
//-------> Calculating MSD for whole plane

void Create_dit_MSD(float *d_data, size_t nTimesamples, size_t nDMs, std::vector<MSD_Data> *dit_MSD, std::vector<MSD_Data> *dit_MSD_BLN, int max_DIT_value, const char *filename, float sigma_constant){
	GpuTimer timer, total_timer;
	double total_time=0, dit_time=0, MSD_time=0, MSD_BLN_time=0;
	int nRest;
	MSD_Data mdtemp;
	size_t decimated_timesamples;
	int DIT_value;
	float *d_sudy, *d_lichy, *d_MSD;
	float h_MSD[MSD_RESULTS_SIZE];
	char str[200];
	cudaMalloc((void **) &d_lichy, (nTimesamples>>1)*nDMs*sizeof(float));
	cudaMalloc((void **) &d_sudy, (nTimesamples>>2)*nDMs*sizeof(float));
	cudaMalloc((void **) &d_MSD, MSD_RESULTS_SIZE*sizeof(float));
	
	total_timer.Start();
	
	
	//----------------------------------------------------------------------------------------
	DIT_value = 1;
	printf("DiT:%d; nTimesamples:%d; decimated_timesamples:%d\n", (int) DIT_value, (int) nTimesamples, (int) (nTimesamples>>1));
	timer.Start();
	MSD_normal(d_data, d_MSD, nTimesamples, nDMs, 0);
	timer.Stop();	MSD_time += timer.Elapsed();
	cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	dit_MSD->push_back(mdtemp);
	
	timer.Start();
	MSD_outlier_rejection(d_data, d_MSD, nDMs, nTimesamples, 0, sigma_constant);
	timer.Stop();	MSD_BLN_time += timer.Elapsed();
	cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	dit_MSD_BLN->push_back(mdtemp);
	
	timer.Stop();
	total_time +=timer.Elapsed();
	//----------------------------------------------------------------------------------------
	
	//----------------------------------------------------------------------------------------
	timer.Start();
	DIT_value = DIT_value*2;
	printf("DiT:%d; nTimesamples:%d; decimated_timesamples:%d\n", (int) DIT_value, (int) nTimesamples, (int) (nTimesamples>>1));
	nRest = GPU_DiT_v2_wrapper(d_data, d_lichy, nDMs, nTimesamples);
	decimated_timesamples = (nTimesamples>>1);
	timer.Stop();	dit_time += timer.Elapsed();
	
	timer.Start();
	MSD_normal(d_lichy, d_MSD, decimated_timesamples, nDMs, nRest);
	timer.Stop();	MSD_time += timer.Elapsed();
	cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	dit_MSD->push_back(mdtemp);
	
	timer.Start();
	MSD_outlier_rejection(d_lichy, d_MSD, nDMs, decimated_timesamples, nRest, sigma_constant);
	timer.Stop();	MSD_BLN_time += timer.Elapsed();
	cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	dit_MSD_BLN->push_back(mdtemp);	
	timer.Stop();
	total_time +=timer.Elapsed();
	//----------------------------------------------------------------------------------------
	
	checkCudaErrors(cudaGetLastError());
	
	for(size_t f=2; f<=max_DIT_value; f++){
		timer.Start();
		DIT_value = DIT_value*2;
		printf("DiT:%d; nTimesamples:%d; decimated_timesamples:%d\n", (int) DIT_value, (int) decimated_timesamples, (int) (decimated_timesamples>>1));
		sprintf(str,"%s_%d", filename, DIT_value);
		if(f%2==0){
			timer.Start();
			nRest = GPU_DiT_v2_wrapper(d_lichy, d_sudy, nDMs, decimated_timesamples);
			timer.Stop();	dit_time += timer.Elapsed();
			if(nRest<0) break;
			decimated_timesamples = (decimated_timesamples>>1);
			
			timer.Start();
			MSD_normal(d_sudy, d_MSD, decimated_timesamples, nDMs, nRest);
			timer.Stop();	MSD_time += timer.Elapsed();
			cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
			mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
			dit_MSD->push_back(mdtemp);
			
			timer.Start();
			MSD_outlier_rejection(d_sudy, d_MSD, nDMs, decimated_timesamples, nRest, sigma_constant);
			timer.Stop();	MSD_BLN_time += timer.Elapsed();
			cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
			mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
			dit_MSD_BLN->push_back(mdtemp);	
		}
		else {
			timer.Start();
			nRest = GPU_DiT_v2_wrapper(d_sudy, d_lichy, nDMs, decimated_timesamples);
			timer.Stop();	dit_time += timer.Elapsed();
			if(nRest<0) break;
			decimated_timesamples = (decimated_timesamples>>1);
			
			timer.Start();
			MSD_normal(d_lichy, d_MSD, decimated_timesamples, nDMs, nRest);
			timer.Stop();	MSD_time += timer.Elapsed();
			cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
			mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
			dit_MSD->push_back(mdtemp);
			
			timer.Start();
			MSD_outlier_rejection(d_lichy, d_MSD, nDMs, decimated_timesamples, nRest, sigma_constant);
			timer.Stop();	MSD_BLN_time += timer.Elapsed();
			cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
			mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
			dit_MSD_BLN->push_back(mdtemp);	
		}
		checkCudaErrors(cudaGetLastError());
	}
	
	total_timer.Stop();
	total_time = total_timer.Elapsed();
	
	printf("Total time: %f; DiT time: %f; MSD time: %f; MSD BLN time: %f;\n", total_time, dit_time, MSD_time, MSD_BLN_time);
	
	checkCudaErrors(cudaFree(d_sudy));
	checkCudaErrors(cudaFree(d_lichy));
	checkCudaErrors(cudaFree(d_MSD));
}

void Create_boxcar_MSD(float *d_data, size_t nTimesamples, size_t nDMs, std::vector<MSD_Data> *boxcar_MSD, std::vector<MSD_Data> *boxcar_MSD_BLN, int max_nTaps, float sigma_constant){
	GpuTimer timer;
	double total_time = 0;
	int nRest;
	MSD_Data mdtemp;
	float *d_boxcar, *d_MSD;
	float h_MSD[MSD_RESULTS_SIZE];
	cudaMalloc((void **) &d_boxcar, nTimesamples*nDMs*sizeof(float));
	cudaMalloc((void **) &d_MSD, MSD_RESULTS_SIZE*sizeof(float));
	
	timer.Start();
	
	MSD_normal(d_data, d_MSD, nDMs, nTimesamples, 0);
	cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	mdtemp.taps = 1; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	boxcar_MSD->push_back(mdtemp);
	
	MSD_outlier_rejection(d_data, d_MSD, nDMs, nTimesamples, 0, sigma_constant);
	cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	mdtemp.taps = 1; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	boxcar_MSD_BLN->push_back(mdtemp);
	
	timer.Stop();
	total_time = total_time + timer.Elapsed();
	printf("DIT value: %d; took %f ms; Total time %fms\n", 1, timer.Elapsed(), total_time);
	
	for(size_t f=2; f<=max_nTaps; f++){
		if( (nTimesamples-f+1)>0 ) {
			timer.Start();
			
			nRest = PD_FIR(d_data, d_boxcar, f, nDMs, nTimesamples);
			
			MSD_normal(d_boxcar, d_MSD, nDMs, nTimesamples, nRest);
			cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
			mdtemp.taps = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
			boxcar_MSD->push_back(mdtemp);
			
			MSD_outlier_rejection(d_boxcar, d_MSD, nDMs, nTimesamples, nRest, sigma_constant);
			cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
			mdtemp.taps = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
			boxcar_MSD_BLN->push_back(mdtemp);
			
			timer.Stop();
			total_time = total_time + timer.Elapsed();
			printf("DIT value: %d; took %f ms; Total time %fms\n", (int) f, timer.Elapsed(), total_time);
		}
	}
	
	checkCudaErrors(cudaGetLastError());
	
	for(size_t f=130; f<=256; f+=4){
		printf("nTimesamples: %d; f: %d; %d\n", nTimesamples, f, nTimesamples-f+1);
		int itemp = (int) (nTimesamples-f+1);
		if( itemp>0 ) {
			timer.Start();
			
			nRest=PPF_L1(d_data, d_boxcar, nDMs, nTimesamples, f);
			
			if(nRest>0){
				MSD_normal(d_boxcar, d_MSD, nDMs, nTimesamples, nRest);
				cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
				mdtemp.taps = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
				boxcar_MSD->push_back(mdtemp);
				
				MSD_outlier_rejection(d_boxcar, d_MSD, nDMs, nTimesamples, nRest, sigma_constant);
				cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
				mdtemp.taps = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
				boxcar_MSD_BLN->push_back(mdtemp);
			}

			timer.Stop();
			total_time = total_time + timer.Elapsed();
			printf("DIT value: %d; took %f ms; Total time %fms\n", (int) f, timer.Elapsed(), total_time);
		}
		checkCudaErrors(cudaGetLastError());
	}
	
	checkCudaErrors(cudaGetLastError());
	
	for(size_t f=272; f<=512; f+=16){
		printf("nTimesamples: %d; f: %d; %d\n", nTimesamples, f, nTimesamples-f+1);
		int itemp = (int) (nTimesamples-f+1);
		if( itemp>0 ) {
			timer.Start();
			
			nRest=PPF_L1(d_data, d_boxcar, nDMs, nTimesamples, f);
			
			if(nRest>0){
				MSD_normal(d_boxcar, d_MSD, nDMs, nTimesamples, nRest);
				cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
				mdtemp.taps = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
				boxcar_MSD->push_back(mdtemp);
				
				MSD_outlier_rejection(d_boxcar, d_MSD, nDMs, nTimesamples, nRest, sigma_constant);
				cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
				mdtemp.taps = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
				boxcar_MSD_BLN->push_back(mdtemp);
			}
			
			timer.Stop();
			total_time = total_time + timer.Elapsed();
			printf("DIT value: %d; took %f ms; Total time %fms\n", (int) f, timer.Elapsed(), total_time);
		}
		checkCudaErrors(cudaGetLastError());
	}
	
	checkCudaErrors(cudaGetLastError());
	
	for(size_t f=544; f<=1024; f+=32){
		printf("nTimesamples: %d; f: %d; %d\n", nTimesamples, f, nTimesamples-f+1);
		int itemp = (int) (nTimesamples-f+1);
		if( itemp>0 ) {
			timer.Start();
			
			nRest=PPF_L1(d_data, d_boxcar, nDMs, nTimesamples, f);

			if(nRest>0){
				MSD_normal(d_boxcar, d_MSD, nDMs, nTimesamples, nRest);
				cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
				mdtemp.taps = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
				boxcar_MSD->push_back(mdtemp);
				
				MSD_outlier_rejection(d_boxcar, d_MSD, nDMs, nTimesamples, nRest, sigma_constant);
				cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
				mdtemp.taps = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
				boxcar_MSD_BLN->push_back(mdtemp);
			}
			
			timer.Stop();
			total_time = total_time + timer.Elapsed();
			printf("DIT value: %d; took %f ms; Total time %fms\n", (int) f, timer.Elapsed(), total_time);
		}
		checkCudaErrors(cudaGetLastError());
	}
	
	checkCudaErrors(cudaGetLastError());

	for(size_t f=1088; f<=2048; f+=64){
		printf("nTimesamples: %d; f: %d; %d\n", nTimesamples, f, nTimesamples-f+1);
		int itemp = (int) (nTimesamples-f+1);
		if( itemp>0 ) {
			timer.Start();
			
			nRest=PPF_L1(d_data, d_boxcar, nDMs, nTimesamples, f);

			if(nRest>0){
				MSD_normal(d_boxcar, d_MSD, nDMs, nTimesamples, nRest);
				cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
				mdtemp.taps = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
				boxcar_MSD->push_back(mdtemp);
				
				MSD_outlier_rejection(d_boxcar, d_MSD, nDMs, nTimesamples, nRest, sigma_constant);
				cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
				mdtemp.taps = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
				boxcar_MSD_BLN->push_back(mdtemp);
			}
			
			timer.Stop();
			total_time = total_time + timer.Elapsed();
			printf("DIT value: %d; took %f ms; Total time %fms\n", (int) f, timer.Elapsed(), total_time);
		}
		checkCudaErrors(cudaGetLastError());
	}
	
	checkCudaErrors(cudaGetLastError());

	for(size_t f=2176; f<=4096; f+=128){
		printf("nTimesamples: %d; f: %d; %d\n", nTimesamples, f, nTimesamples-f+1);
		int itemp = (int) (nTimesamples-f+1);
		if( itemp>0 ) {
			timer.Start();
			
			nRest=PPF_L1(d_data, d_boxcar, nDMs, nTimesamples, f);

			if(nRest>0){		
				MSD_normal(d_boxcar, d_MSD, nDMs, nTimesamples, nRest);
				cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
				mdtemp.taps = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
				boxcar_MSD->push_back(mdtemp);
				
				MSD_outlier_rejection(d_boxcar, d_MSD, nDMs, nTimesamples, nRest, sigma_constant);
				cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
				mdtemp.taps = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
				boxcar_MSD_BLN->push_back(mdtemp);
			}
			
			timer.Stop();
			total_time = total_time + timer.Elapsed();
			printf("DIT value: %d; took %f ms; Total time %fms\n", (int) f, timer.Elapsed(), total_time);
		}
		checkCudaErrors(cudaGetLastError());
	}
	
	checkCudaErrors(cudaGetLastError());
	
	for(size_t f=4352; f<=8192; f+=256){
		printf("nTimesamples: %d; f: %d; %d\n", nTimesamples, f, nTimesamples-f+1);
		int itemp = (int) (nTimesamples-f+1);
		if( itemp>0 ) {
			timer.Start();
			
			nRest=PPF_L1(d_data, d_boxcar, nDMs, nTimesamples, f);
			
			if(nRest>0){
				MSD_normal(d_boxcar, d_MSD, nDMs, nTimesamples, nRest);
				cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
				mdtemp.taps = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
				boxcar_MSD->push_back(mdtemp);
				
				MSD_outlier_rejection(d_boxcar, d_MSD, nDMs, nTimesamples, nRest, sigma_constant);
				cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
				mdtemp.taps = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
				boxcar_MSD_BLN->push_back(mdtemp);
			}
			
			timer.Stop();
			total_time = total_time + timer.Elapsed();
			printf("DIT value: %d; took %f ms; Total time %fms\n", (int) f, timer.Elapsed(), total_time);
		}
		checkCudaErrors(cudaGetLastError());
	}
	
	checkCudaErrors(cudaGetLastError());
	checkCudaErrors(cudaFree(d_boxcar));
	checkCudaErrors(cudaFree(d_MSD));
}

void Export_MSD_data(std::vector<MSD_Data> h_dit_MSD, std::vector<MSD_Data> h_dit_MSD_BLN, std::vector<MSD_Data> h_boxcar_MSD, std::vector<MSD_Data> h_boxcar_MSD_BLN, char *filename){
	std::ofstream FILEOUT;
	FILEOUT.open (filename, std::ofstream::out);

	for(size_t f=0; f<h_dit_MSD.size(); f++){
		FILEOUT << (int) h_dit_MSD[f].taps << " " << h_dit_MSD[f].mean << " " << h_dit_MSD[f].sd << " " << "1" << std::endl;
	}
	FILEOUT << std::endl;
	FILEOUT << std::endl;
	for(size_t f=0; f<h_dit_MSD_BLN.size(); f++){
		FILEOUT << (int) h_dit_MSD_BLN[f].taps << " " << h_dit_MSD_BLN[f].mean << " " << h_dit_MSD_BLN[f].sd << " " << "2" << std::endl;
	}
	FILEOUT << std::endl;
	FILEOUT << std::endl;
	
	for(size_t f=0; f<h_boxcar_MSD.size(); f++){
		FILEOUT << (int) h_boxcar_MSD[f].taps << " " << h_boxcar_MSD[f].mean << " " << h_boxcar_MSD[f].sd << " " << "3" << std::endl;
	}
	FILEOUT << std::endl;
	FILEOUT << std::endl;
	for(size_t f=0; f<h_boxcar_MSD_BLN.size(); f++){
		FILEOUT << (int) h_boxcar_MSD_BLN[f].taps << " " << h_boxcar_MSD_BLN[f].mean << " " << h_boxcar_MSD_BLN[f].sd << " " << "4" << std::endl;
	}
	FILEOUT << std::endl;
	FILEOUT << std::endl;	
	
	FILEOUT.close();
}

void Calculate_MSD_data(float *output_buffer, size_t nTimesamples, size_t nDMs, float sigma_constant, int inBin, float dm_low, float dm_high, float tstart){
	char filename[200];
	int max_DIT_value = 13;
	int max_nTaps = 128;
	std::vector<MSD_Data> h_dit_MSD;
	std::vector<MSD_Data> h_dit_MSD_BLN;
	std::vector<MSD_Data> h_boxcar_MSD;
	std::vector<MSD_Data> h_boxcar_MSD_BLN;
	
	size_t free_mem,total_mem;
	cudaMemGetInfo(&free_mem,&total_mem);
	printf("Memory available: %f; output_buffer size: %f;\n", (double) free_mem/(1024.0*1024.0), ((double) nDMs*nTimesamples*sizeof(float))/(1024.0*1024.0));
	
	Create_dit_MSD(output_buffer, nTimesamples/inBin, nDMs, &h_dit_MSD, &h_dit_MSD_BLN, max_DIT_value, filename, sigma_constant);
	Create_boxcar_MSD(output_buffer, nTimesamples/inBin, nDMs, &h_boxcar_MSD, &h_boxcar_MSD_BLN, max_nTaps, sigma_constant);
	
	
	sprintf(filename,"MSD_test-t_%.2f-dm_%.2f-%.2f.dat", tstart, dm_low, dm_high);
	Export_MSD_data(h_dit_MSD, h_dit_MSD_BLN, h_boxcar_MSD, h_boxcar_MSD_BLN, filename);
}

//-------> Calculating MSD for whole plane
//---------------------------------------------------------------------------------



// Extend this to arbitrary size plans
void Create_PD_plan(std::vector<PulseDetection_plan> *PD_plan, std::vector<int> *BC_widths, int nDMs, int nTimesamples){
	int Elements_per_block, itemp, nRest;
	PulseDetection_plan PDmp;
	
	if(BC_widths->size()>0){
		PDmp.shift        = 0;
		PDmp.output_shift = 0;
		PDmp.startTaps    = 0;
		PDmp.iteration    = 0;
		
		PDmp.decimated_timesamples = nTimesamples;
		PDmp.dtm = (nTimesamples>>(PDmp.iteration+1));
		PDmp.dtm = PDmp.dtm - (PDmp.dtm&1);
		
		PDmp.nBoxcars = BC_widths->operator[](0);
		Elements_per_block = PD_NTHREADS*2 - PDmp.nBoxcars;
		itemp = PDmp.decimated_timesamples;
		PDmp.nBlocks = itemp/Elements_per_block;
		nRest = itemp - PDmp.nBlocks*Elements_per_block;
		if(nRest>0) PDmp.nBlocks++;
		PDmp.unprocessed_samples = PDmp.nBoxcars + 6;
		if(PDmp.decimated_timesamples<PDmp.unprocessed_samples) PDmp.nBlocks=0;
		PDmp.total_ut = PDmp.unprocessed_samples;
		
		
		PD_plan->push_back(PDmp);
		
		for(int f=1; f< (int) BC_widths->size(); f++){
			// These are based on previous values of PDmp
			PDmp.shift        = PDmp.nBoxcars/2;
			PDmp.output_shift = PDmp.output_shift + PDmp.decimated_timesamples;
			PDmp.startTaps    = PDmp.startTaps + PDmp.nBoxcars*(1<<PDmp.iteration);
			PDmp.iteration    = PDmp.iteration + 1;
			
			// Definition of new PDmp values
			PDmp.decimated_timesamples = PDmp.dtm;
			PDmp.dtm = (nTimesamples>>(PDmp.iteration+1));
			PDmp.dtm = PDmp.dtm - (PDmp.dtm&1);
			
			PDmp.nBoxcars = BC_widths->operator[](f);
			Elements_per_block=PD_NTHREADS*2 - PDmp.nBoxcars;
			itemp = PDmp.decimated_timesamples;
			PDmp.nBlocks = itemp/Elements_per_block;
			nRest = itemp - PDmp.nBlocks*Elements_per_block;
			if(nRest>0) PDmp.nBlocks++;
			PDmp.unprocessed_samples = PDmp.unprocessed_samples/2 + PDmp.nBoxcars + 6; //
			if(PDmp.decimated_timesamples<PDmp.unprocessed_samples) PDmp.nBlocks=0;
			PDmp.total_ut = PDmp.unprocessed_samples*(1<<PDmp.iteration);
			
			PD_plan->push_back(PDmp);
		}
	}
}


int Get_max_iteration(int max_boxcar_width, std::vector<int> *BC_widths, int *max_width_performed){
	int startTaps, iteration;
	
	startTaps = 0;
	iteration = 0;
	for(int f=0; f<(int) BC_widths->size(); f++){
		startTaps = startTaps + BC_widths->operator[](f)*(1<<f);
		if(startTaps>=max_boxcar_width) {
			iteration = f+1;
			break;
		}
	}
	
	if(max_boxcar_width>startTaps) {
		iteration=(int) BC_widths->size();
	}
	
	*max_width_performed=startTaps;
	return(iteration);
}


void analysis_GPU(float *h_peak_list, size_t *peak_pos, size_t max_peak_size, int i, float tstart, int t_processed, int inBin, int outBin, int *maxshift, int max_ndms, int *ndms, float cutoff, float sigma_constant, float max_boxcar_width_in_sec, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp, int candidate_algorithm, int enable_sps_baselinenoise){
	int max_boxcar_width = (int) (max_boxcar_width_in_sec/tsamp);
	int max_width_performed=0;
	//unsigned long int j;
	unsigned long int vals;
	int nTimesamples = t_processed;
	int nDMs = ndms[i];
	int temp_peak_pos;
	
	//----------------------------------------------
	//--- MSD profile of the data
	//Calculate_MSD_data(output_buffer, nTimesamples, nDMs, sigma_constant, inBin, dm_low[i], dm_high[i], tstart);
	//---------------------------------------------<
	

	// Calculate the total number of values
	vals = (unsigned long int) ( nDMs*nTimesamples );
	

	double total_time=0, partial_time=0, MSD_time=0, SPDT_time=0, PF_time=0;
	
	
	//float max, min, threshold;
	int max_iteration;
	int t_BC_widths[10]={PD_MAXTAPS,16,16,16,8,8,8,8,8,8};
	std::vector<int> BC_widths(t_BC_widths,t_BC_widths+sizeof(t_BC_widths)/sizeof(int));
	std::vector<PulseDetection_plan> PD_plan;

	//---------------------------------------------------------------------------
	//----------> GPU part
	printf("\n----------> GPU analysis part\n");
	printf("     Dimensions nDMs:%d; nTimesamples:%d; inBin:%d; outBin:%d; maxshift:%d; \n", ndms[i], t_processed, inBin, outBin, *maxshift);
	GpuTimer total_timer, timer;
	total_timer.Start();
	
	//float h_MSD[3];
	float *d_MSD;
	checkCudaErrors(cudaGetLastError());
	if ( cudaSuccess != cudaMalloc((void**) &d_MSD, sizeof(float)*3)) {printf("Allocation error!\n"); exit(201);}
	
	/*
	//-------------- CPU check
	float *h_temp, *h_boxcar;
	double signal_mean, signal_sd;
	h_temp = (float *)malloc( ((size_t) nDMs*nTimesamples)*sizeof(float));
	h_boxcar = (float *)malloc( ((size_t) nDMs*nTimesamples)*sizeof(float));
	memset(h_temp, 0.0, ((size_t) nDMs*nTimesamples)*sizeof(float));
	memset(h_boxcar, 0.0, ((size_t) nDMs*nTimesamples)*sizeof(float));
	cudaMemcpy( h_temp, output_buffer, ((size_t) nDMs*nTimesamples)*sizeof(float), cudaMemcpyDeviceToHost);
	MSD_Kahan(h_temp, nDMs, nTimesamples, 0, &signal_mean, &signal_sd);
	printf("MSD_kahan: after 1 tap   Mean: %e, Standard deviation: %e;\n",signal_mean, signal_sd);
	Calculate_FIR(h_temp, h_boxcar, 2, nDMs, nTimesamples, 0);
	MSD_Kahan(h_boxcar, nDMs, nTimesamples, 0, &signal_mean, &signal_sd);
	printf("MSD_kahan: after 1 tap   Mean: %e, Standard deviation: %e;\n",signal_mean, signal_sd);	
	free(h_temp);
	free(h_boxcar);
	//-------------- CPU check
	*/
	
	/*
	//-------------- One Call linear approximation
	timer.Start(); 
	MSD_linear_approximation(output_buffer, d_MSD, PD_MAXTAPS, nDMs, nTimesamples, 0);
	timer.Stop();
	partial_time = timer.Elapsed(); 
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	printf("     MSD linear approximation: Mean: %f, Stddev: %f, modifier: %f\n", h_MSD[0], h_MSD[1], h_MSD[2]);
	#ifdef GPU_ANALYSIS_DEBUG
	printf("     One kernel took:%f ms\n", partial_time); 
	#endif
	//-------------- One Call linear approximation
	*/
	
	/*
	//-------------- Base level noise point-wise
	timer.Start(); 
	MSD_outlier_rejection(output_buffer, d_MSD, nDMs, nTimesamples, 0, sigma_constant);
	timer.Stop();
	partial_time = timer.Elapsed(); 
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	printf("     MSD BLN point-wise: Mean: %f, Stddev: %f, modifier: %f\n", h_MSD[0], h_MSD[1], h_MSD[2]);
	#ifdef GPU_ANALYSIS_DEBUG
	printf("     MSD BLN point-wise kernel took:%f ms\n", partial_time); 
	#endif
	//-------------- Base level noise point-wise
	*/

	/*
	//-------------- BLN_LA
	timer.Start(); 
	MSD_BLN_LA_pw_normal(output_buffer, d_MSD, nDMs, nTimesamples, PD_MAXTAPS, 0, sigma_constant);
	timer.Stop();
	partial_time = timer.Elapsed(); 
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	printf("     MSD BLN linear approximation: Mean: %f, Stddev: %f, modifier: %f\n", h_MSD[0], h_MSD[1], h_MSD[2]);
	#ifdef GPU_ANALYSIS_DEBUG
	printf("     BLN LA took:%f ms\n", partial_time); 
	#endif
	//-------------- BLN_LA
	*/
	
	/*
	//-------------- Base level noise grid
	timer.Start(); 
	MSD_BLN_grid(output_buffer, d_MSD, 32, 32, nDMs, nTimesamples, 0, sigma_constant);
	timer.Stop();
	partial_time = timer.Elapsed(); 
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	printf("     MSD BLN grid: Mean: %f, Stddev: %f, modifier: %f\n", h_MSD[0], h_MSD[1], h_MSD[2]);
	#ifdef GPU_ANALYSIS_DEBUG
	printf("     MSD BLN grid kernel took:%f ms\n", partial_time); 
	#endif
	//-------------- Base level noise grid
	*/
	
	size_t free_mem,total_mem;
	cudaMemGetInfo(&free_mem,&total_mem);
	printf("     Memory required by boxcar filters:%0.3f MB\n",(4.5*vals*sizeof(float) + 2*vals*sizeof(ushort))/(1024.0*1024) );
	printf("     Memory available:%0.3f MB \n", ((float) free_mem)/(1024.0*1024.0) );
	
	std::vector<int> DM_list;
	unsigned long int max_timesamples=(free_mem*0.95)/(5.5*sizeof(float) + 2*sizeof(ushort));
	int DMs_per_cycle = max_timesamples/nTimesamples;
	int nRepeats, nRest, DM_shift, itemp, local_max_list_size;//BC_shift,
	
	itemp = (int) (DMs_per_cycle/THR_WARPS_PER_BLOCK);
	DMs_per_cycle = itemp*THR_WARPS_PER_BLOCK;
	
	nRepeats = nDMs/DMs_per_cycle;
	nRest = nDMs - nRepeats*DMs_per_cycle;
	local_max_list_size = (DMs_per_cycle*nTimesamples)/4;
	
	for(int f=0; f<nRepeats; f++) DM_list.push_back(DMs_per_cycle);
	if(nRest>0) DM_list.push_back(nRest);
	
	printf("     SPS will run %d batches each containing %d DM trials. Remainder %d DM trials\n", (int) DM_list.size(), DMs_per_cycle, nRest);
	
	
	max_iteration = Get_max_iteration(max_boxcar_width/inBin, &BC_widths, &max_width_performed);
	printf("     Selected iteration:%d; maximum boxcar width requested:%d; maximum boxcar width performed:%d;\n", max_iteration, max_boxcar_width/inBin, max_width_performed);
	Create_PD_plan(&PD_plan, &BC_widths, 1, nTimesamples);
	
	if(DM_list.size()>0){
		DMs_per_cycle = DM_list[0];
		
		float *d_peak_list;
		if ( cudaSuccess != cudaMalloc((void**) &d_peak_list, sizeof(float)*DMs_per_cycle*nTimesamples)) printf("Allocation error! peaks\n");
		
		float *d_decimated;
		if ( cudaSuccess != cudaMalloc((void **) &d_decimated,  sizeof(float)*(((DMs_per_cycle*nTimesamples)/2)+PD_MAXTAPS) )) printf("Allocation error! dedispered\n");
		
		float *d_boxcar_values;
		if ( cudaSuccess != cudaMalloc((void **) &d_boxcar_values,  sizeof(float)*DMs_per_cycle*nTimesamples)) printf("Allocation error! boxcars\n");
		
		float *d_output_SNR;
		if ( cudaSuccess != cudaMalloc((void **) &d_output_SNR, sizeof(float)*2*DMs_per_cycle*nTimesamples)) printf("Allocation error! SNR\n");
		
		ushort *d_output_taps;
		if ( cudaSuccess != cudaMalloc((void **) &d_output_taps, sizeof(ushort)*2*DMs_per_cycle*nTimesamples)) printf("Allocation error! taps\n");
		
		int *gmem_peak_pos;
		cudaMalloc((void**) &gmem_peak_pos, 1*sizeof(int));
		cudaMemset((void*) gmem_peak_pos, 0, sizeof(int));
		
		DM_shift = 0;
		for(int f=0; f<DM_list.size(); f++) {
			//-------------- SPS BLN
			timer.Start();
			//PD_SEARCH_LONG_BLN(&output_buffer[DM_shift*nTimesamples], d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, &PD_plan, max_iteration, DM_list[f], nTimesamples);
			//PD_SEARCH_LONG_BLN_EACH(&output_buffer[DM_shift*nTimesamples], d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, &PD_plan, max_iteration, DM_list[f], nTimesamples, sigma_constant);
			//PD_SEARCH_LONG_LINAPPROX(&output_buffer[DM_shift*nTimesamples], d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, &PD_plan, max_iteration, DM_list[f], nTimesamples);
			if(enable_sps_baselinenoise){
				PD_SEARCH_LONG_BLN_LINAPPROX_EACH(&output_buffer[DM_shift*nTimesamples], d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, &PD_plan, max_iteration, DM_list[f], nTimesamples, sigma_constant, &MSD_time, &SPDT_time);
			}
			else {
				PD_SEARCH_LONG_LINAPPROX_EACH(&output_buffer[DM_shift*nTimesamples], d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, &PD_plan, max_iteration, DM_list[f], nTimesamples, &MSD_time, &SPDT_time);
			}
			//
			timer.Stop();
			partial_time = timer.Elapsed();
			#ifdef GPU_ANALYSIS_DEBUG
			printf("PD_SEARCH took:%f ms\n", partial_time);
			#endif
			//-------------- SPS BLN
			
			checkCudaErrors(cudaGetLastError());
			
			#ifdef GPU_ANALYSIS_DEBUG
			printf("BC_shift:%d; DMs_per_cycle:%d; f*DMs_per_cycle:%d; max_iteration:%d;\n", DM_shift*nTimesamples, DM_list[f], DM_shift, max_iteration);
			#endif
			
			if(candidate_algorithm==1){
				//-------------- Thresholding
				timer.Start();
				THRESHOLD(d_output_SNR, d_output_taps, d_peak_list, gmem_peak_pos, cutoff, DM_list[f], nTimesamples, DM_shift, &PD_plan, max_iteration, local_max_list_size);
				timer.Stop();
				partial_time = timer.Elapsed();
				PF_time = timer.Elapsed();
				#ifdef GPU_ANALYSIS_DEBUG
				printf("THR_WARP took:%f ms\n", partial_time);
				#endif
				//-------------- Thresholding
			}
			else {
				//-------------- Peak finding
				timer.Start();
				PEAK_FIND(d_output_SNR, d_output_taps, d_peak_list, DM_list[f], nTimesamples, cutoff, local_max_list_size, gmem_peak_pos, DM_shift, &PD_plan, max_iteration);
				timer.Stop();
				partial_time = timer.Elapsed();
				PF_time = timer.Elapsed();
				#ifdef GPU_ANALYSIS_DEBUG
				printf("PEAK_FIND took:%f ms\n", partial_time);
				#endif
				//-------------- Peak finding
			}
			
			checkCudaErrors(cudaGetLastError());
			
			checkCudaErrors(cudaMemcpy(&temp_peak_pos, gmem_peak_pos, sizeof(int), cudaMemcpyDeviceToHost));
			#ifdef GPU_ANALYSIS_DEBUG
			printf("temp_peak_pos:%d; host_pos:%zu; max:%zu; local_max:%d;\n", temp_peak_pos, (*peak_pos), max_peak_size, local_max_list_size);
			#endif
			if( temp_peak_pos>=local_max_list_size ) {
				printf("     Maximum list size reached! Increase list size or increase sigma cutoff.\n");
				temp_peak_pos=local_max_list_size;
			}
			if( ((*peak_pos) + temp_peak_pos)<max_peak_size){
				checkCudaErrors(cudaMemcpy(&h_peak_list[(*peak_pos)*4], d_peak_list, temp_peak_pos*4*sizeof(float), cudaMemcpyDeviceToHost));
				*peak_pos = (*peak_pos) + temp_peak_pos;
			}
			else printf("Error peak list is too small!\n");
			

			//---------> Old thresholding code.
			//#ifdef OLD_THRESHOLD
			//#endif
			//---------> Old thresholding code.

			DM_shift = DM_shift + DM_list[f];
			cudaMemset((void*) gmem_peak_pos, 0, sizeof(int));
		}
		
		//------------------------> Output
		#pragma omp parallel for
		for (int count = 0; count < (*peak_pos); count++){
			h_peak_list[4*count]     = h_peak_list[4*count]*dm_step[i] + dm_low[i];
			h_peak_list[4*count + 1] = h_peak_list[4*count + 1]*tsamp + tstart;
		}
        
		FILE *fp_out;
		char filename[200];
		
		if(candidate_algorithm==1){
			if((*peak_pos)>0){
				sprintf(filename, "analysed-t_%.2f-dm_%.2f-%.2f.dat", tstart, dm_low[i], dm_high[i]);
				if (( fp_out = fopen(filename, "wb") ) == NULL)	{
					fprintf(stderr, "Error opening output file!\n");
					exit(0);
				}
				fwrite(h_peak_list, (*peak_pos)*sizeof(float), 4, fp_out);
				fclose(fp_out);
			}
		}
		else {
			if((*peak_pos)>0){
				sprintf(filename, "peak_analysed-t_%.2f-dm_%.2f-%.2f.dat", tstart, dm_low[i], dm_high[i]);
				if (( fp_out = fopen(filename, "wb") ) == NULL)	{
					fprintf(stderr, "Error opening output file!\n");
					exit(0);
				}
				fwrite(h_peak_list, (*peak_pos)*sizeof(float), 4, fp_out);
				fclose(fp_out);
			}
		}
		//------------------------> Output
		
		cudaFree(d_peak_list);
		cudaFree(d_boxcar_values);
		cudaFree(d_decimated);
		cudaFree(d_output_SNR);
		cudaFree(d_output_taps);
		cudaFree(gmem_peak_pos);

	}
	else printf("Error not enough memory to search for pulses\n");

	total_timer.Stop();
	total_time = total_timer.Elapsed();
	printf("\n     TOTAL TIME OF SPS:%f ms\n", total_time);
	printf("\n     MSD_time: %f ms; SPDT time: %f ms; Peak find time: %f ms;\n", MSD_time, SPDT_time, PF_time);
	printf("----------<\n\n");
	total_time += partial_time;

	cudaFree(d_MSD);
	//----------> GPU part
	//---------------------------------------------------------------------------
	
}
