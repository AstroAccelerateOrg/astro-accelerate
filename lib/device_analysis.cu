//#define GPU_ANALYSIS_DEBUG
//#define MSD_PLANE_EXPORT
//#define GPU_PARTIAL_TIMER
#define GPU_TIMER


#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "headers/params.h"

#include "headers/device_BC_plan.h"
#include "headers/device_peak_find.h"
#include "headers/device_MSD_Configuration.h"
#include "headers/device_MSD.h"
#include "headers/device_SPS_long.h"
#include "headers/device_threshold.h"
#include "headers/device_single_FIR.h"

#include "timer.h"

//TODO:
// Make BC_plan for arbitrary long pulses, by reusing last element in the plane



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

void MSD_on_GPU(float *h_input, float *d_input, float *d_MSD, float *signal_mean, float *signal_sd, float *signal_mean_bln, float *signal_sd_bln, float *signal_mean_bl_bln, float *signal_sd_bl_bln, int nDMs, int nTimesamples, int offset, float OR_sigma_multiplier, float *MSD_limited_time, float *MSD_BLN_pw_time, float *MSD_BLN_grid_time){
	GpuTimer timer;
	float h_MSD[3];
	cudaMemcpy( d_input, h_input, ((size_t) nDMs*nTimesamples)*sizeof(float), cudaMemcpyHostToDevice);
	
	timer.Start();
	MSD_normal(d_MSD, d_input, nTimesamples, nDMs, offset);
	timer.Stop();
	(*MSD_limited_time) = timer.Elapsed();
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	(*signal_mean) = h_MSD[0];
	(*signal_sd)   = h_MSD[1];
	
	
	timer.Start();
	MSD_outlier_rejection(d_MSD, d_input, nTimesamples, nDMs, offset, OR_sigma_multiplier);
	timer.Stop();
	(*MSD_BLN_pw_time) = timer.Elapsed();
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	(*signal_mean_bln) = h_MSD[0];
	(*signal_sd_bln)   = h_MSD[1];
	
	
	timer.Start();
	MSD_grid_outlier_rejection(d_input, d_MSD, 32, 32, nTimesamples, nDMs, offset, OR_sigma_multiplier);
	timer.Stop();
	(*MSD_BLN_grid_time) = timer.Elapsed();
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost); 
	(*signal_mean_bl_bln) = h_MSD[0];
	(*signal_sd_bl_bln)   = h_MSD[1];
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
//-------> Calculating MSD of the whole plane for different boxcar widths

//TODO: Write a kernel which would calculate MSD for DIT=1 and DIT=2 and outputs plane for DIT=4. This will save some memory.
//		Provide a way of choosing type of interpolation.
//		Throw out the garbage code.

void do_MSD_both(float *d_MSD, float *d_MSD_OR, float *d_input, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, double *MSD_time, double *MSD_OR_time){
	GpuTimer timer;
	
	timer.Start();
	MSD_normal(d_MSD, d_input, nTimesamples, nDMs, offset);
	timer.Stop();	(*MSD_time) += timer.Elapsed();
	
	timer.Start();
	MSD_outlier_rejection(d_MSD_OR, d_input, nTimesamples, nDMs, offset, OR_sigma_multiplier);
	timer.Stop();	(*MSD_OR_time) += timer.Elapsed();
}


void Do_MSD(float *d_MSD, float *d_input, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, int MSD_type){
	if(MSD_type){
		MSD_outlier_rejection(d_MSD, d_input, nTimesamples, nDMs, offset, OR_sigma_multiplier);
	}
	else {
		MSD_normal(d_MSD, d_input, nTimesamples, nDMs, offset);
	}
}


void MSD_of_input_plane(float *d_MSD_DIT, float *d_data, std::vector<int> *h_MSD_DIT_widths, size_t nTimesamples, size_t nDMs, int nDecimations, int max_width_performed, float OR_sigma_multiplier, int MSD_type){
	GpuTimer timer, total_timer;
	double total_time=0, dit_time=0, MSD_time=0;
	int nRest;
	size_t decimated_timesamples;
	int DIT_value;
	float *d_sudy, *d_lichy;
	cudaMalloc((void **) &d_lichy, (nTimesamples>>1)*nDMs*sizeof(float));
	cudaMalloc((void **) &d_sudy, (nTimesamples>>2)*nDMs*sizeof(float));
	
	total_timer.Start();
	//----------------------------------------------------------------------------------------
	//-------- DIT = 1
	DIT_value = 1;
	
	timer.Start();
	Do_MSD(d_MSD_DIT, d_data, nTimesamples, nDMs, 0, OR_sigma_multiplier, MSD_type);
	timer.Stop();	MSD_time += timer.Elapsed();
	h_MSD_DIT_widths->push_back(DIT_value);

	#ifdef GPU_ANALYSIS_DEBUG
		float h_MSD[MSD_RESULTS_SIZE];
		printf("    MSD format: [ mean ; StDev ; nElements ]\n");
		checkCudaErrors(cudaMemcpy(h_MSD, d_MSD_DIT, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost));
		printf("    DiT:%d; nTimesamples:%d; decimated_timesamples:%d; MSD:[%f; %f; %f]\n", (int) DIT_value, (int) nTimesamples, (int) (nTimesamples>>1), h_MSD[0], h_MSD[1], h_MSD[2]);
	#endif
	//----------------------------------------------------------------------------------------
	
	checkCudaErrors(cudaGetLastError());
	
	//----------------------------------------------------------------------------------------
	//-------- DIT = 2
	timer.Start();
	DIT_value = DIT_value*2;
	
	nRest = GPU_DiT_v2_wrapper(d_data, d_lichy, nDMs, nTimesamples);
	decimated_timesamples = (nTimesamples>>1);
	timer.Stop();	dit_time += timer.Elapsed();

	timer.Start();
	Do_MSD(&d_MSD_DIT[MSD_RESULTS_SIZE], d_lichy, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, MSD_type);
	timer.Stop();	MSD_time += timer.Elapsed();
	h_MSD_DIT_widths->push_back(DIT_value);
	
	#ifdef GPU_ANALYSIS_DEBUG
		cudaMemcpy(h_MSD, &d_MSD_DIT[MSD_RESULTS_SIZE], MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
		printf("    DiT:%d; nTimesamples:%d; decimated_timesamples:%d; MSD:[%f; %f; %f]\n", (int) DIT_value, (int) nTimesamples, (int) (nTimesamples>>1), h_MSD[0], h_MSD[1], h_MSD[2]);
	#endif
	//----------------------------------------------------------------------------------------
	
	checkCudaErrors(cudaGetLastError());
	
	//----------------------------------------------------------------------------------------
	//-------- DIT > 2
	for(size_t f=2; f<=nDecimations; f++){
		timer.Start();
		DIT_value = DIT_value*2;
		if(f%2==0){
			timer.Start();
			nRest = GPU_DiT_v2_wrapper(d_lichy, d_sudy, nDMs, decimated_timesamples);
			timer.Stop();	dit_time += timer.Elapsed();
			if(nRest<0) break;
			decimated_timesamples = (decimated_timesamples>>1);

			timer.Start();
			Do_MSD(&d_MSD_DIT[f*MSD_RESULTS_SIZE], d_sudy, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, MSD_type);
			timer.Stop();	MSD_time += timer.Elapsed();
		}
		else {
			timer.Start();
			nRest = GPU_DiT_v2_wrapper(d_sudy, d_lichy, nDMs, decimated_timesamples);
			timer.Stop();	dit_time += timer.Elapsed();
			if(nRest<0) break;
			decimated_timesamples = (decimated_timesamples>>1);

			timer.Start();
			Do_MSD(&d_MSD_DIT[f*MSD_RESULTS_SIZE], d_lichy, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, MSD_type);
			timer.Stop();	MSD_time += timer.Elapsed();
		}
		h_MSD_DIT_widths->push_back(DIT_value);
		
		#ifdef GPU_ANALYSIS_DEBUG
			cudaMemcpy(h_MSD, &d_MSD_DIT[f*MSD_RESULTS_SIZE], MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
			printf("    DiT:%d; nTimesamples:%d; decimated_timesamples:%d; MSD:[%f; %f; %f]\n", (int) DIT_value, (int) decimated_timesamples, (int) (decimated_timesamples>>1), h_MSD[0], h_MSD[1], h_MSD[2]);
		#endif
		
		checkCudaErrors(cudaGetLastError());
	}
	//----------------------------------------------------------------------------------------
	
	checkCudaErrors(cudaGetLastError());
	
	//----------------------------------------------------------------------------------------
	//-------- Boxcar for last boxcar width if needed
	/*
	if(DIT_value<max_width_performed){
		DIT_value = (DIT_value>>1);
		decimated_timesamples = nTimesamples/DIT_value;
		int nTaps = max_width_performed/DIT_value;
		if(max_width_performed%DIT_value!=0) nTaps++;
		
		if(nDecimations%2==0){
			nRest = PPF_L1(d_lichy, d_sudy, nDMs, decimated_timesamples, nTaps);

			checkCudaErrors(cudaGetLastError());
			
			timer.Start();
			Do_MSD(&d_MSD_DIT[(nDecimations+1)*MSD_RESULTS_SIZE], d_sudy, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, MSD_type);
			timer.Stop();	MSD_time += timer.Elapsed();
		}
		else {
			nRest = PPF_L1(d_sudy, d_lichy, nDMs, decimated_timesamples, nTaps);
			
			checkCudaErrors(cudaGetLastError());
			
			timer.Start();
			Do_MSD(&d_MSD_DIT[(nDecimations+1)*MSD_RESULTS_SIZE], d_lichy, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, MSD_type);
			timer.Stop();	MSD_time += timer.Elapsed();
		}
		h_MSD_DIT_widths->push_back(DIT_value*nTaps);

		#ifdef GPU_ANALYSIS_DEBUG
			printf("    Performing additional boxcar: nTaps: %d; max_width_performed: %d; DIT_value/2: %d;\n", nTaps, max_width_performed, DIT_value);
			checkCudaErrors(cudaMemcpy(h_MSD, &d_MSD_DIT[(nDecimations+1)*MSD_RESULTS_SIZE], MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost));
			printf("    DIT: %d; MSD:[%f; %f; %f]\n", DIT_value*nTaps, h_MSD[0], h_MSD[1], h_MSD[2]);
		#endif		
	}
	*/
	//----------------------------------------------------------------------------------------
	
	checkCudaErrors(cudaGetLastError());
	
	checkCudaErrors(cudaFree(d_sudy));
	checkCudaErrors(cudaFree(d_lichy));
	
	total_timer.Stop();
	total_time = total_timer.Elapsed();
	
	#ifdef GPU_PARTIAL_TIMER
		printf("    MSD of input plane: Total time: %f ms; DiT time: %f ms; MSD time: %f ms;\n", total_time, dit_time, MSD_time);
	#endif
}


void Interpolate_linear(float *mean, float *StDev, float desired_width, float *h_MSD_DIT, std::vector<int> *h_MSD_DIT_widths){
	int MSD_DIT_size = h_MSD_DIT_widths->size();
	int position = (int) floorf(log2f((float) desired_width));
	
	float width1 = h_MSD_DIT_widths->operator[](position);
	float mean1 = h_MSD_DIT[(position)*MSD_RESULTS_SIZE];
	float StDev1 = h_MSD_DIT[(position)*MSD_RESULTS_SIZE +1];
	
	if(position == MSD_DIT_size-1 && width1==(int) desired_width) {
		(*mean) = mean1;
		(*StDev) = StDev1;
	}
	else {
		float width2 = h_MSD_DIT_widths->operator[](position+1);
		float distance_in_width = width2 - width1;
		
		float mean2 = h_MSD_DIT[(position+1)*MSD_RESULTS_SIZE];
		float distance_in_mean = mean2 - mean1;
		
		float StDev2 = h_MSD_DIT[(position+1)*MSD_RESULTS_SIZE +1];
		float distance_in_StDev = StDev2 - StDev1;
	
		(*mean) = mean1 + (distance_in_mean/distance_in_width)*((float) desired_width - width1);
		(*StDev) = StDev1 + (distance_in_StDev/distance_in_width)*((float) desired_width - width1);
	}
}


void Interpolate_square(float *mean, float *StDev, float desired_width, float *h_MSD_DIT, std::vector<int> *h_MSD_DIT_widths){
	int MSD_DIT_size = h_MSD_DIT_widths->size();
	int position = (int) floorf(log2f((float) desired_width));
	
	if(position == MSD_DIT_size-2) position--;
	if(position == MSD_DIT_size-1 && h_MSD_DIT_widths->operator[](position)==(int) desired_width) {
		(*mean)  = h_MSD_DIT[(position)*MSD_RESULTS_SIZE];
		(*StDev) = h_MSD_DIT[(position)*MSD_RESULTS_SIZE +1];
	}
	else {
		float w = desired_width;
		
		float w0 = h_MSD_DIT_widths->operator[](position);
		float mean0  = h_MSD_DIT[(position)*MSD_RESULTS_SIZE];
		float StDev0 = h_MSD_DIT[(position)*MSD_RESULTS_SIZE +1];
		
		float w1 = h_MSD_DIT_widths->operator[](position+1);
		float mean1  = h_MSD_DIT[(position+1)*MSD_RESULTS_SIZE];
		float StDev1 = h_MSD_DIT[(position+1)*MSD_RESULTS_SIZE +1];
		
		float w2 = h_MSD_DIT_widths->operator[](position+2);
		float mean2  = h_MSD_DIT[(position+2)*MSD_RESULTS_SIZE];
		float StDev2 = h_MSD_DIT[(position+2)*MSD_RESULTS_SIZE +1];
		
		float a0 = ((w - w1)*(w - w2))/((w0 - w1)*(w0 - w2));
		float a1 = ((w - w0)*(w - w2))/((w1 - w0)*(w1 - w2));
		float a2 = ((w - w0)*(w - w1))/((w2 - w0)*(w2 - w1));
		
		(*mean)  = a0*mean0 + a1*mean1 + a2*mean2;
		(*StDev) = a0*StDev0 + a1*StDev1 + a2*StDev2;
	}
}


void Export_MSD_plane(const char *filename, float *h_MSD_DIT, std::vector<int> *h_MSD_DIT_widths, float *h_MSD_interpolated, std::vector<int> *h_boxcar_widths, int max_width_performed) {
	char str[200];
	std::ofstream FILEOUT;
	int MSD_INTER_SIZE = 2;
	
	sprintf(str,"%s_DIT.dat", filename);
	FILEOUT.open (str, std::ofstream::out);
	for(size_t f=0; f<(int) h_MSD_DIT_widths->size(); f++){
		FILEOUT << (int) h_MSD_DIT_widths->operator[](f) << " " << h_MSD_DIT[f*MSD_RESULTS_SIZE] << " " << h_MSD_DIT[f*MSD_RESULTS_SIZE + 1] << std::endl;
	}
	FILEOUT.close();
	
	sprintf(str,"%s_Interpolated.dat", filename);
	FILEOUT.open (str, std::ofstream::out);
	for(size_t f=0; f<(int) h_boxcar_widths->size(); f++){
		if(h_boxcar_widths->operator[](f)<=max_width_performed)
			FILEOUT << (int) h_boxcar_widths->operator[](f) << " " << h_MSD_interpolated[f*MSD_INTER_SIZE] << " " << h_MSD_interpolated[f*MSD_INTER_SIZE + 1] << std::endl;
	}
	FILEOUT.close();
}


void Interpolate_MSD_values(float *d_MSD_interpolated, float *d_MSD_DIT, std::vector<int> *h_MSD_DIT_widths, int nMSDs, std::vector<int> *h_boxcar_widths, int max_width_performed, const char *filename){
	#ifdef GPU_PARTIAL_TIMER
	GpuTimer timer;
	timer.Start();
	#endif
	
	int MSD_INTER_SIZE = 2;
	float *h_MSD_DIT, *h_MSD_interpolated;
	int nWidths = (int) h_boxcar_widths->size();
	h_MSD_DIT = new float[nMSDs*MSD_RESULTS_SIZE];
	h_MSD_interpolated = new float[nWidths*MSD_INTER_SIZE];
	
	checkCudaErrors(cudaMemcpy(h_MSD_DIT, d_MSD_DIT, nMSDs*MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost));
	
	for(int f=0; f<nWidths; f++){
		if(h_boxcar_widths->operator[](f)<=max_width_performed) {
			float mean, StDev;
			Interpolate_square(&mean, &StDev, (float) h_boxcar_widths->operator[](f), h_MSD_DIT, h_MSD_DIT_widths);
			h_MSD_interpolated[f*MSD_INTER_SIZE] = mean;
			h_MSD_interpolated[f*MSD_INTER_SIZE+1] = StDev;
		}
	}
	
	#ifdef MSD_PLANE_EXPORT
		Export_MSD_plane(filename, h_MSD_DIT, h_MSD_DIT_widths, h_MSD_interpolated, h_boxcar_widths, max_width_performed);
	#endif
	
	checkCudaErrors(cudaMemcpy(d_MSD_interpolated, h_MSD_interpolated, nWidths*MSD_INTER_SIZE*sizeof(float), cudaMemcpyHostToDevice));
	
	delete[] h_MSD_DIT;
	delete[] h_MSD_interpolated;
	
	#ifdef GPU_PARTIAL_TIMER
	timer.Stop();
	printf("    Interpolation step took %f ms;\n", timer.Elapsed());
	#endif
}

//-------------------------------------------------------------------------<











//---------------------------------------------------------------------------------
//-------> Calculating MSD for whole plane via DIT and boxcar

void Create_dit_MSD(float *d_data, size_t nTimesamples, size_t nDMs, std::vector<MSD_Data> *dit_MSD, std::vector<MSD_Data> *dit_MSD_BLN, int max_DIT_value, const char *filename, float OR_sigma_multiplier){
	GpuTimer timer, total_timer;
	double total_time=0, dit_time=0, MSD_time=0, MSD_OR_time=0;
	int nRest;
	MSD_Data mdtemp;
	size_t decimated_timesamples;
	int DIT_value;
	float *d_sudy, *d_lichy, *d_MSD, *d_MSD_OR;
	float h_MSD[MSD_RESULTS_SIZE];
	char str[200];
	checkCudaErrors(cudaMalloc((void **) &d_lichy, (nTimesamples>>1)*nDMs*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **) &d_sudy, (nTimesamples>>2)*nDMs*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **) &d_MSD, MSD_RESULTS_SIZE*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **) &d_MSD_OR, MSD_RESULTS_SIZE*sizeof(float)));
	
	total_timer.Start();
	
	checkCudaErrors(cudaGetLastError());
	//----------------------------------------------------------------------------------------
	DIT_value = 1;
	printf("DiT:%d; nTimesamples:%d; decimated_timesamples:%d\n", (int) DIT_value, (int) nTimesamples, (int) (nTimesamples>>1));
	
	do_MSD_both(d_MSD, d_MSD_OR, d_data, nTimesamples, nDMs, 0, OR_sigma_multiplier, &MSD_time, &MSD_OR_time);
	
	checkCudaErrors(cudaGetLastError());
	
	checkCudaErrors(cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost));
	mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	dit_MSD->push_back(mdtemp);
	
	checkCudaErrors(cudaMemcpy(h_MSD, d_MSD_OR, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost));
	mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	dit_MSD_BLN->push_back(mdtemp);
	//----------------------------------------------------------------------------------------
	
	//----------------------------------------------------------------------------------------
	timer.Start();
	DIT_value = DIT_value*2;
	printf("DiT:%d; nTimesamples:%d; decimated_timesamples:%d\n", (int) DIT_value, (int) nTimesamples, (int) (nTimesamples>>1));
	nRest = GPU_DiT_v2_wrapper(d_data, d_lichy, nDMs, nTimesamples);
	decimated_timesamples = (nTimesamples>>1);
	timer.Stop();	dit_time += timer.Elapsed();

	do_MSD_both(d_MSD, d_MSD_OR, d_lichy, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, &MSD_time, &MSD_OR_time);
	
	cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	dit_MSD->push_back(mdtemp);
	
	cudaMemcpy(h_MSD, d_MSD_OR, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	dit_MSD_BLN->push_back(mdtemp);
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

			do_MSD_both(d_MSD, d_MSD_OR, d_sudy, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, &MSD_time, &MSD_OR_time);
			
			cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
			mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
			dit_MSD->push_back(mdtemp);
			
			cudaMemcpy(h_MSD, d_MSD_OR, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
			mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
			dit_MSD_BLN->push_back(mdtemp);
		}
		else {
			timer.Start();
			nRest = GPU_DiT_v2_wrapper(d_sudy, d_lichy, nDMs, decimated_timesamples);
			timer.Stop();	dit_time += timer.Elapsed();
			if(nRest<0) break;
			decimated_timesamples = (decimated_timesamples>>1);

			do_MSD_both(d_MSD, d_MSD_OR, d_lichy, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, &MSD_time, &MSD_OR_time);
			
			cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
			mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
			dit_MSD->push_back(mdtemp);
			
			cudaMemcpy(h_MSD, d_MSD_OR, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
			mdtemp.taps = DIT_value; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
			dit_MSD_BLN->push_back(mdtemp);
		}
		checkCudaErrors(cudaGetLastError());
	}
	
	total_timer.Stop();
	total_time = total_timer.Elapsed();
	
	printf("Total time: %f; DiT time: %f; MSD time: %f; MSD BLN time: %f;\n", total_time, dit_time, MSD_time, MSD_OR_time);
	
	checkCudaErrors(cudaFree(d_sudy));
	checkCudaErrors(cudaFree(d_lichy));
	checkCudaErrors(cudaFree(d_MSD));
	checkCudaErrors(cudaFree(d_MSD_OR));
}

void Create_boxcar_MSD(float *d_data, size_t nTimesamples, size_t nDMs, std::vector<MSD_Data> *boxcar_MSD, std::vector<MSD_Data> *boxcar_MSD_BLN, int max_nTaps, float OR_sigma_multiplier){
	/*
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
	
	MSD_outlier_rejection(d_data, d_MSD, nDMs, nTimesamples, 0, OR_sigma_multiplier);
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
			
			MSD_outlier_rejection(d_boxcar, d_MSD, nDMs, nTimesamples, nRest, OR_sigma_multiplier);
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
				
				MSD_outlier_rejection(d_boxcar, d_MSD, nDMs, nTimesamples, nRest, OR_sigma_multiplier);
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
				
				MSD_outlier_rejection(d_boxcar, d_MSD, nDMs, nTimesamples, nRest, OR_sigma_multiplier);
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
				
				MSD_outlier_rejection(d_boxcar, d_MSD, nDMs, nTimesamples, nRest, OR_sigma_multiplier);
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
				
				MSD_outlier_rejection(d_boxcar, d_MSD, nDMs, nTimesamples, nRest, OR_sigma_multiplier);
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
				
				MSD_outlier_rejection(d_boxcar, d_MSD, nDMs, nTimesamples, nRest, OR_sigma_multiplier);
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
				
				MSD_outlier_rejection(d_boxcar, d_MSD, nDMs, nTimesamples, nRest, OR_sigma_multiplier);
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
	*/
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

void Calculate_MSD_data(float *output_buffer, size_t nTimesamples, size_t nDMs, float OR_sigma_multiplier, int inBin, float dm_low, float dm_high, float tstart){
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
	
	Create_dit_MSD(output_buffer, nTimesamples/inBin, nDMs, &h_dit_MSD, &h_dit_MSD_BLN, max_DIT_value, filename, OR_sigma_multiplier);
	Create_boxcar_MSD(output_buffer, nTimesamples/inBin, nDMs, &h_boxcar_MSD, &h_boxcar_MSD_BLN, max_nTaps, OR_sigma_multiplier);
	
	
	sprintf(filename,"MSD_test-t_%.2f-dm_%.2f-%.2f.dat", tstart, dm_low, dm_high);
	Export_MSD_data(h_dit_MSD, h_dit_MSD_BLN, h_boxcar_MSD, h_boxcar_MSD_BLN, filename);
}

//-----------------------------------------------------------------------------------<




void Create_list_of_boxcar_widths(std::vector<int> *boxcar_widths, std::vector<int> *BC_widths){
	int DIT_value, DIT_factor, width;
	DIT_value = 1;
	DIT_factor = 2;
	width = 0;
	for(int f=0; f<(int) BC_widths->size(); f++){
		for(int b=0; b<BC_widths->operator[](f); b++){
			width = width + DIT_value;
			boxcar_widths->push_back(width);
		}
		DIT_value = DIT_value*DIT_factor;
	}
}

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


void analysis_GPU(float *h_peak_list, size_t *peak_pos, size_t max_peak_size, int i, float tstart, int t_processed, int inBin, int outBin, int *maxshift, int max_ndms, int *ndms, float cutoff, float OR_sigma_multiplier, float max_boxcar_width_in_sec, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp, int candidate_algorithm, int enable_sps_baselinenoise){
	//--------> Task
	int max_boxcar_width = (int) (max_boxcar_width_in_sec/tsamp);
	int max_width_performed=0;
	size_t vals;
	int nTimesamples = t_processed;
	int nDMs = ndms[i];
	int temp_peak_pos;
	
	//--------> Benchmarking
	double total_time=0, MSD_time=0, SPDT_time=0, PF_time=0;
	
	//--------> Other
	size_t free_mem,total_mem;
	char filename[200];
	
	//----------------------------------------------
	//--- 
	cudaMemGetInfo(&free_mem,&total_mem);
	//printf("Memory available: %f; output_buffer size: %f;\n", (double) free_mem/(1024.0*1024.0), ((double) nDMs*nTimesamples*sizeof(float))/(1024.0*1024.0));
	//Calculate_MSD_data(output_buffer, nTimesamples, nDMs, OR_sigma_multiplier, inBin, dm_low[i], dm_high[i], tstart);
	//---------------------------------------------<
	
	
	// Calculate the total number of values
	vals = ((size_t) nDMs)*((size_t) nTimesamples);
	
	
	//float max, min, threshold;
	int max_iteration;
	int t_BC_widths[10]={PD_MAXTAPS,16,16,16,8,8,8,8,8,8};
	std::vector<int> BC_widths(t_BC_widths,t_BC_widths+sizeof(t_BC_widths)/sizeof(int));
	std::vector<PulseDetection_plan> PD_plan;

	//---------------------------------------------------------------------------
	//----------> GPU part
	printf("\n----------> GPU analysis part\n");
	printf("  Dimensions nDMs:%d; nTimesamples:%d; inBin:%d; outBin:%d; maxshift:%d; \n", ndms[i], t_processed, inBin, outBin, *maxshift);
	GpuTimer total_timer, timer;
	total_timer.Start();
	
	
	float *d_MSD;
	checkCudaErrors(cudaGetLastError());
	if ( cudaSuccess != cudaMalloc((void**) &d_MSD, sizeof(float)*3)) {printf("Allocation error!\n"); exit(201);}
	
	
	cudaMemGetInfo(&free_mem,&total_mem);
	printf("  Memory required by boxcar filters:%0.3f MB\n",(4.5*vals*sizeof(float) + 2*vals*sizeof(ushort))/(1024.0*1024) );
	printf("  Memory available:%0.3f MB \n", ((float) free_mem)/(1024.0*1024.0) );
	
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
	
	if( (int) DM_list.size() > 1 ) 
		printf("  SPS will run %d batches each containing %d DM trials. Remainder %d DM trials\n", (int) DM_list.size(), DMs_per_cycle, nRest);
	else 
		printf("  SPS will run %d batch containing %d DM trials.\n", (int) DM_list.size(), nRest);
	
	
	max_iteration = Get_max_iteration(max_boxcar_width/inBin, &BC_widths, &max_width_performed);
	printf("  Selected iteration:%d; maximum boxcar width requested:%d; maximum boxcar width performed:%d;\n", max_iteration, max_boxcar_width/inBin, max_width_performed);
	Create_PD_plan(&PD_plan, &BC_widths, 1, nTimesamples);
	std::vector<int> h_boxcar_widths;
	Create_list_of_boxcar_widths(&h_boxcar_widths, &BC_widths);
	
	
	/*
	//---------------------------------------------
	//------- Calculation of MSD for whole plane
	timer.Start();
	
	int MSD_INTER_SIZE = 2;
	float *d_MSD_DIT, *d_MSD_interpolated;
	int nDIT_widths, nboxcar_Widths, nDecimations;
	std::vector<int> h_MSD_DIT_widths;
	
	nDecimations = ((int) floorf(log2f((float)max_width_performed))) + 1;
	nDIT_widths = nDecimations + 1;
	
	
	nboxcar_Widths = h_boxcar_widths.size();
	
	cudaMalloc((void **) &d_MSD_DIT, nDIT_widths*MSD_RESULTS_SIZE*sizeof(float));
	cudaMalloc((void **) &d_MSD_interpolated, nboxcar_Widths*MSD_INTER_SIZE*sizeof(float));
	
	MSD_of_input_plane(d_MSD_DIT, output_buffer, &h_MSD_DIT_widths, nTimesamples, nDMs, nDecimations, max_width_performed, OR_sigma_multiplier, enable_sps_baselinenoise);
	
	#ifdef GPU_ANALYSIS_DEBUG
		printf("    Number of calculated MSD values: %d; number of interpolated MSD values: %d;\n",nDIT_widths, nboxcar_Widths);
	#endif
	
	sprintf(filename,"MSD_interpolate_test-t_%.2f-dm_%.2f-%.2f", tstart, dm_low[i], dm_high[i]);
	Interpolate_MSD_values(d_MSD_interpolated, d_MSD_DIT, &h_MSD_DIT_widths, nDIT_widths, &h_boxcar_widths, max_width_performed, filename);
	
	timer.Stop();
	MSD_time = timer.Elapsed();
	#ifdef GPU_PARTIAL_TIMER
	printf("    MSD_plane time: %f ms\n", MSD_time);
	#endif
	//printf("Konec\n");
	//------- Calculation of MSD for whole plane
	//---------------------------------------------
	*/
	
	
	
	
	//-------------------------------------------------------------------------
	//------------ Using MSD_plane_profile
	size_t MSD_profile_size_in_bytes, MSD_DIT_profile_size_in_bytes, workarea_size_in_bytes;
	cudaMemGetInfo(&free_mem,&total_mem);
	Get_MSD_plane_profile_memory_requirements(&MSD_profile_size_in_bytes, &MSD_DIT_profile_size_in_bytes, &workarea_size_in_bytes, nTimesamples, nDMs, &h_boxcar_widths);
	double dit_time, MSD_only_time;
	float *d_MSD_interpolated;
	float *d_MSD_DIT = NULL;
	float *temporary_workarea;
	cudaMalloc((void **) &d_MSD_interpolated, MSD_profile_size_in_bytes);
	cudaMalloc((void **) &temporary_workarea, workarea_size_in_bytes);
	
	MSD_plane_profile(d_MSD_interpolated, output_buffer, d_MSD_DIT, temporary_workarea, false, nTimesamples, nDMs, &h_boxcar_widths, tstart, dm_low[i], dm_high[i], OR_sigma_multiplier, enable_sps_baselinenoise, false, &total_time, &dit_time, &MSD_time);
	
	#ifdef GPU_PARTIAL_TIMER
		printf("    MSD time: Total: %f ms; DIT: %f ms; MSD: %f ms;\n", MSD_time, dit_time, MSD_only_time);
	#endif
	
	cudaFree(temporary_workarea);
	//------------ Using MSD_plane_profile
	//-------------------------------------------------------------------------	
	
	
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
			//-------------- SPDT
			timer.Start();
			SPDT_search_long_MSD_plane(&output_buffer[DM_shift*nTimesamples], d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD_interpolated, &PD_plan, max_iteration, nTimesamples, DM_list[f]);
			timer.Stop();
			SPDT_time += timer.Elapsed();
			#ifdef GPU_PARTIAL_TIMER
			printf("    SPDT took:%f ms\n", timer.Elapsed());
			#endif
			//-------------- SPDT
			
			checkCudaErrors(cudaGetLastError());
			
			#ifdef GPU_ANALYSIS_DEBUG
			printf("    BC_shift:%d; DMs_per_cycle:%d; f*DMs_per_cycle:%d; max_iteration:%d;\n", DM_shift*nTimesamples, DM_list[f], DM_shift, max_iteration);
			#endif
			
			if(candidate_algorithm==1){
				//-------------- Thresholding
				timer.Start();
				THRESHOLD(d_output_SNR, d_output_taps, d_peak_list, gmem_peak_pos, cutoff, DM_list[f], nTimesamples, DM_shift, &PD_plan, max_iteration, local_max_list_size);
				timer.Stop();
				PF_time += timer.Elapsed();
				#ifdef GPU_PARTIAL_TIMER
				printf("    Thresholding took:%f ms\n", timer.Elapsed());
				#endif
				//-------------- Thresholding
			}
			else {
				//-------------- Peak finding
				timer.Start();
				PEAK_FIND(d_output_SNR, d_output_taps, d_peak_list, DM_list[f], nTimesamples, cutoff, local_max_list_size, gmem_peak_pos, DM_shift, &PD_plan, max_iteration);
				timer.Stop();
				PF_time = timer.Elapsed();
				#ifdef GPU_PARTIAL_TIMER
				printf("    Peak finding took:%f ms\n", timer.Elapsed());
				#endif
				//-------------- Peak finding
			}
			
			checkCudaErrors(cudaGetLastError());
			
			checkCudaErrors(cudaMemcpy(&temp_peak_pos, gmem_peak_pos, sizeof(int), cudaMemcpyDeviceToHost));
			#ifdef GPU_ANALYSIS_DEBUG
			printf("    temp_peak_pos:%d; host_pos:%zu; max:%zu; local_max:%d;\n", temp_peak_pos, (*peak_pos), max_peak_size, local_max_list_size);
			#endif
			if( temp_peak_pos>=local_max_list_size ) {
				printf("    Maximum list size reached! Increase list size or increase sigma cutoff.\n");
				temp_peak_pos=local_max_list_size;
			}
			if( ((*peak_pos) + temp_peak_pos)<max_peak_size){
				checkCudaErrors(cudaMemcpy(&h_peak_list[(*peak_pos)*4], d_peak_list, temp_peak_pos*4*sizeof(float), cudaMemcpyDeviceToHost));
				*peak_pos = (*peak_pos) + temp_peak_pos;
			}
			else printf("Error peak list is too small!\n");

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
		cudaFree(d_MSD_DIT);
		cudaFree(d_MSD_interpolated);

	}
	else printf("Error not enough memory to search for pulses\n");

	total_timer.Stop();
	total_time = total_timer.Elapsed();
	#ifdef GPU_TIMER
	printf("\n  TOTAL TIME OF SPS:%f ms\n", total_time);
	printf("  MSD_time: %f ms; SPDT time: %f ms; Candidate selection time: %f ms;\n", MSD_time, SPDT_time, PF_time);
	printf("----------<\n\n");
	#endif

	cudaFree(d_MSD);
	//----------> GPU part
	//---------------------------------------------------------------------------
	
}
