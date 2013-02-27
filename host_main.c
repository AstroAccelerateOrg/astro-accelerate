#include "headers_mains.h"
#include <xmmintrin.h>

#include "device_dedisperse.h"
#include "device_init.h"
#include "device_malloc.h"
#include "device_load_data.h"
#include "device_save_data.h"
#include "device_free_memory.h"

#include "host_bin.h"
#include "host_analysis.h"
#include "host_init.h"
#include "host_help.h"

int main(int argc, char *argv[])
{

	//{{{ Declare and initialise variables
	
	int c, i;
	int kernel_type = 0;

	// Telescope parameters
	int nchans, nsamp, nbits, nsamples, nifs;

	double start, end, tstart, tsamp, fch1, foff;

	// Get date from cli and open the binary data file
	FILE *fp, *fp_out;
	
	// De-dispersion parameters
	int tdms;

	float dm_step;
	float dm_low;
	
	//}}}
	
	// Get de-dispersion parameters
	
	init(argc, argv, &fp, &tdms, &dm_step, &dm_low);


	//{{{ Read in the significant bits of the header
	
	char *string = (char *) malloc(80 * sizeof(char));
	int nchar;
	int nbytes = sizeof(int);

	while(1) {
	
		strcpy(string,"ERROR");
		fread(&nchar, sizeof(int), 1, fp);

		if (feof(fp)) exit(0);

		if (nchar > 1 && nchar < 80) {
			
			fread(string, nchar, 1, fp);
			string[nchar] = '\0';
			// For debugging only
			//printf("\n%d\t%s", nchar, string), fflush(stdout);
			nbytes += nchar;

			if (strcmp(string,"HEADER_END") == 0) break;
		
			if (strcmp(string, "tsamp")           == 0) {
				fread(&tsamp,  sizeof(double), 1, fp);
			} else if (strcmp(string, "tstart")   == 0) {
				fread(&tstart, sizeof(double), 1, fp);
			} else if (strcmp(string, "fch1")     == 0) {
				fread(&fch1,   sizeof(double), 1, fp);
			} else if (strcmp(string, "foff")     == 0) {
				fread(&foff,   sizeof(double), 1, fp);
			} else if (strcmp(string, "nchans")   == 0) {
				fread(&nchans,    sizeof(int), 1, fp);
			} else if (strcmp(string, "nifs")     == 0) {
				fread(&nifs,      sizeof(int), 1, fp);
			} else if (strcmp(string, "nbits")    == 0) {
				fread(&nbits,     sizeof(int), 1, fp);
			} else if (strcmp(string, "nsamples") == 0) {
				fread(&nsamples,  sizeof(int), 1, fp);
			}
		}
	}

	// For Debugging only
	
	printf("\ntsamp:\t\t%lf", tsamp);
	printf("\ntstart:\t\t%lf", tstart);
	printf("\nfch1:\t\t%lf", fch1);
	printf("\nfoff:\t\t%lf", foff);
	printf("\nnchans:\t\t%d", nchans);
	printf("\nnifs:\t\t%d", nifs);
	printf("\nnbits:\t\t%d", nbits);

	// Check that we are working with one IF channel
	if ( nifs != 1) {
		printf("\nERROR!! Can only work with one IF channel!\n");
		return (1);
	}

	// Note the where the header ends
	fpos_t file_loc;
	fgetpos(fp, &file_loc);

	//}}}

	//{{{ Load in the raw data from the input file and transpose

	printf("\n\n!!LOAD DATA!!");

	start = omp_get_wtime();
	
	// Allocate a tempory buffer to store a line of frequency data
	//float *temp_buffer = (float *) _mm_malloc(nchans * sizeof(float),ALIGNSIZE);
	float *temp_buffer = (float *) malloc(nchans * sizeof(float));

	// Count how many time samples we have
	unsigned long int total_data = 0;
	while(!feof(fp)) {
		fread(temp_buffer, (nbits/8), nchans, fp);
		total_data++;
	}
	nsamp = total_data;
		
	printf("\n\nnsamp:\t\t%d", nsamp);

	// Move the file pointer back to the end of the header
	fsetpos(fp, &file_loc);

	// Allocate memory for input array
	
	float *input_buffer = NULL;
	size_t inputsize =  nsamp * nchans * sizeof(float);
	input_buffer = (float *) _mm_malloc(inputsize,ALIGNSIZE);
	//cudaHostAlloc((void**)&input_buffer,inputsize,cudaHostAllocDefault);

	printf("\n\nInput size:\t%d MB", (int) (inputsize / 1024 / 1024));

	// Read in the data, transpose it and store it in the input buffer

	total_data = 0;
	while(!feof(fp)) {
		if(fread(temp_buffer, (nbits/8), nchans, fp) != nchans) break;
		for(c = 0; c < nchans; c++) {
			input_buffer[c * nsamp + total_data] = temp_buffer[c];
		}
		total_data++;
	}
	
	end = omp_get_wtime();
	printf("\n\nLoad omp time:\t%.16g s\n", end - start);
	
	//}}}

	//{{{ Calculate time independent dm shifts
	
	float *dmshifts = (float *) malloc(nchans * sizeof(float));
	for (c = 0; c < nchans; c++) {
		dmshifts[c] = 4148.741601 * ((1.0 / (fch1 + (foff * c)) / (fch1 + (foff * c))) - (1.0 / fch1 / fch1));
	}
	
	//}}}
	
	//{{{ Initialise the GPU 
	
	printf("\n\n!!INITIALISE GPU!!");

	start = omp_get_wtime();
	
	init_gpu(nchans, dmshifts);

	end = omp_get_wtime();

	printf("\n\ninit_gpu omp time:\t\t\t%.16g s", end - start);	

	//}}}
	
	//{{{ Calculate time binning parameters
	
	// This method relies on defining points when nsamps is a multiple of
	// nchans - bin on the diagonal or a fraction of it.

	float	max_dm		= tdms * dm_step;
	float	fmin		= (fch1 + (foff * nchans));
	float	fmin_squared	= fmin*fmin;
	float	fmax_squared	= fch1*fch1;
	float	bin_dm		= ((nchans * tsamp)/(4148.741601))*((fmin_squared*fmax_squared)/(fmax_squared - fmin_squared));
	int	nbins		= ceil(max_dm/bin_dm);

	printf("\n\nMax DM:\t\t%f", max_dm);
	printf("\nBin Point:\t%f", bin_dm);
	printf("\nNumber of bins:\t%d", nbins);

	//}}}

	//{{{ Check f,t line length
	// Check to see if the threadblock will load a shared memory line that
	// is long enough for the algorithm to run without an out of bounds
	// access...
	
	float	n;
	int	ndms;

	modff((((bin_dm / dm_step) + DIVINDM) / DIVINDM), &n);
	ndms = (int)(n * DIVINDM);

	// Calculate lineshift for shared memory
	int lineshift = (dm_step * ((ndms - 1) * 4148.741601 * ((1.0 / (fch1 + (foff * (nchans - 1))) / (fch1 + (foff * (nchans - 1)))) - (1.0 / fch1 / fch1)) - (ndms - DIVINDM) 
				* 4148.741601 * ((1.0 / (fch1 + (foff * (nchans - DIVINDM))) / (fch1 + (foff * (nchans - DIVINDM)))) - (1.0 / fch1 / fch1))))/tsamp;
	printf("\n\nlineshift:\t%d", lineshift);

	if(((NUMREG - 1)*(DIVINT) + lineshift) > (DIVINT * DIVINDM) && kernel_type == 0) {
		printf("\nERROR!! Your thread block is not large enough to load all of the data needed into shared memory\n\n");
		return(1);
	}

	//}}}

	//{{{ Prepare memory layout
	
	float *output_buffer = NULL;
	float *bin_buffer = NULL;

	float *d_input = NULL;
	float *d_output = NULL;

	//{{{ Calculate maxshift, the number of dms for this bin and
	//the highest value of dm to be calculated in this bin

	modff((((bin_dm / dm_step) + DIVINDM) / DIVINDM), &n);
	ndms = (int)(n * DIVINDM);
	
	//}}}
	
	size_t outputsize = nsamp * ndms * sizeof(float);
	output_buffer = (float *) malloc(outputsize);
	//cudaHostAlloc((void**)&output_buffer, outputsize, cudaHostAllocDefault);

	size_t binsize =  (nsamp * nchans * sizeof(float)) + 1;
	bin_buffer = (float *) malloc(binsize);
	//cudaHostAlloc((void**)&bin_buffer, binsize, cudaHostAllocDefault);

	cudaMalloc((void **) &d_input, inputsize);
	cudaMalloc((void **) &d_output, outputsize);
	cudaMemset(d_output, 0, outputsize);

	//}}}

	//{{{ Start time binning
	
	start = omp_get_wtime();
	//#pragma omp parallel for private(max_dm,dm_low,dm_step,n,ndms,tsamp,nsamp,dmshifts,input_buffer)
	for(i = 1; i <= nbins; i++) {
		
		//{{{ Calculate maxshift, the number of dms for this bin and
		//the highest value of dm to be calculated in this bin

		if(i == nbins) {
			modff(((((max_dm - dm_low)/dm_step) + DIVINDM) / DIVINDM), &n);
			ndms = (int)(n * DIVINDM);
		} else {
			modff((((bin_dm / dm_step) + DIVINDM) / DIVINDM), &n);
			ndms = (int)(n * DIVINDM);
		}
	
		int maxshift = ((dm_low + dm_step * (ndms - 1)) * dmshifts[nchans - 1])/tsamp;
		printf("\nmaxshift:\t%d", maxshift), fflush(stdout);

		float dm_high = dm_low +(ndms*dm_step);
		
		//}}}
	
		load_data(d_input, input_buffer, inputsize, nsamp, maxshift);
		
		printf("\n\n!!DEDISPERSE ON GPU!! Bin: %d DM: %f to %f", i, dm_low, dm_high);
		dedisperse(inputsize, d_input, outputsize, d_output, nchans, nsamp, maxshift, dm_low, ndms, kernel_type, tsamp, dm_step);
		
		size_t binsize =  ((nsamp/powf(2,i)) + 1) * nchans * sizeof(float);
		bin(binsize, bin_buffer, inputsize, input_buffer, nchans, nsamp);
		cudaThreadSynchronize();
		save_data(d_output, output_buffer, outputsize);
		cudaMemset(d_output, 0, outputsize);
		
		//{{{ Perform the analysis for this bin

		analysis(nsamp, maxshift, ndms, outputsize, output_buffer, dm_low, dm_high, dm_step, tsamp);
		
		//}}}

		dm_low = dm_high;
		input_buffer = bin_buffer;
		inputsize = binsize;
		nsamp = nsamp/2;
		tsamp = tsamp*2;
		dm_step = dm_step*2;
		
	}
	end = omp_get_wtime();
	printf("\n\ndedisperse omp time:\t\t\t\t\t%.16g s\n", end - start);

	//}}}
	
	cudaFreeHost(output_buffer);
	free_device_memory(d_input);
	free_device_memory(d_output);

	fclose(fp);

	free(string);
	free(temp_buffer);
	cudaFreeHost(input_buffer);

	return 0;
}
