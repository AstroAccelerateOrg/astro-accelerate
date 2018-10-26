#include <stdio.h>
#include "host_get_file_data.hpp"
/* Note we send in a pointer to the file pointer becuase this function needs to update the position of the file pointer
 */

void get_file_data(FILE **fp, int *nchans, int *nsamples, int *nsamp, int *nifs, int *nbits, float *tsamp, float *tstart, float *fch1, float *foff)
{

	fpos_t file_loc;

	char *string = (char *) malloc(80 * sizeof(char));

	int nchar;
	int nbytes = sizeof(int);

	unsigned long int total_data;

	double temp;

	while (1)
	{

		strcpy(string, "ERROR");
		if (fread(&nchar, sizeof(int), 1, *fp) != 1)
		{
			fprintf(stderr, "\nError while reading file\n");
			exit(0);
		}
		if (feof(*fp))
			exit(0);

		if (nchar > 1 && nchar < 80)
		{
			if (fread(string, nchar, 1, *fp) != 1)
			{
				fprintf(stderr, "\nError while reading file\n");
				exit(0);
			}

			string[nchar] = '\0';
			// For debugging only
			printf("\n%d\t%s", nchar, string), fflush(stdout);
			nbytes += nchar;

			if (strcmp(string, "HEADER_END") == 0)
				break;

			if (strcmp(string, "tsamp") == 0)
			{
				if (fread(&temp, sizeof(double), 1, *fp) != 1)
				{
					fprintf(stderr, "\nError while reading file\n");
					exit(0);
				}
				*tsamp = (float) temp;
			}
			else if (strcmp(string, "tstart") == 0)
			{
				if (fread(&temp, sizeof(double), 1, *fp) != 1)
				{
					fprintf(stderr, "\nError while reading file\n");
					exit(0);
				}
				*tstart = (float) temp;
			}
			else if (strcmp(string, "fch1") == 0)
			{
				if (fread(&temp, sizeof(double), 1, *fp) != 1)
				{
					fprintf(stderr, "\nError while reading file\n");
					exit(0);
				}
				*fch1 = (float) temp;
			}
			else if (strcmp(string, "foff") == 0)
			{
				if (fread(&temp, sizeof(double), 1, *fp) != 1)
				{
					fprintf(stderr, "\nError while reading file\n");
					exit(0);
				}
				*foff = (float) temp;
			}
			else if (strcmp(string, "nchans") == 0)
			{
				if (fread(nchans, sizeof(int), 1, *fp) != 1)
				{
					fprintf(stderr, "\nError while reading file\n");
					exit(0);
				}
			}
			else if (strcmp(string, "nifs") == 0)
			{
				if (fread(nifs, sizeof(int), 1, *fp) != 1)
				{
					fprintf(stderr, "\nError while reading file\n");
					exit(0);
				}
			}
			else if (strcmp(string, "nbits") == 0)
			{
				if (fread(nbits, sizeof(int), 1, *fp) != 1)
				{
					fprintf(stderr, "\nError while reading file\n");
					exit(0);
				}
			}
			else if (strcmp(string, "nsamples") == 0)
			{
				if (fread(nsamples, sizeof(int), 1, *fp) != 1)
				{
					fprintf(stderr, "\nError while reading file\n");
					exit(0);
				}
			}
		}
	}

	// Check that we are working with one IF channel
	if (*nifs != 1)
	{
		printf("\nERROR!! Can only work with one IF channel!\n");
		exit(1);
	}

	// Getting number of time-samples based on file size
	unsigned long int data_start = ftell(*fp);
	if (fseek(*fp, 0, SEEK_END) != 0) {
		printf("\nERROR!! Failed to seek to the end of data file\n");
		exit(1);
	}
	unsigned long int exp_total_data = ftell(*fp);
	exp_total_data = exp_total_data - data_start;
	fseek(*fp, data_start, SEEK_SET);

	if (( *nbits ) == 32) {
		*nsamp = exp_total_data/((*nchans)*4);
	}
	else if (( *nbits ) == 16) {
		*nsamp = exp_total_data/((*nchans)*2);
	}
	else if (( *nbits ) == 8) {
		*nsamp = exp_total_data/((*nchans));
	}
	else if (( *nbits ) == 4) {
		*nsamp = 2*exp_total_data/((*nchans));
	}
	else if (( *nbits ) == 2) {
		*nsamp = 4*exp_total_data/((*nchans));
	}
	else if (( *nbits ) == 1) {
		*nsamp = 8*exp_total_data/((*nchans));
	}
	else {
		printf("\n\n======================= ERROR ==========================\n");
		printf(    " Currently this code only runs with 1, 2, 4 8 and 16 bit data \n");
		printf(  "\n========================================================\n");
		exit(0);
	}

}
