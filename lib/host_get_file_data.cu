#include <stdio.h>

/* Note we send in a pointer to the file pointer becuase this function needs to update the position of the file pointer
 */

void get_file_data(FILE **fp, int *nchans, int *nsamples, int *nsamp, int *nifs, int *nbits, float *tsamp, float *tstart, float *fch1, float *foff)
{

	fpos_t file_loc;

	char *string = (char *) malloc(80 * sizeof(char));

	int nchar;
	int nbytes = sizeof(int);

	long int total_data;

	double temp;

	while (1)
	{

		strcpy(string, "ERROR");
		fread(&nchar, sizeof(int), 1, *fp);

		if (feof(*fp))
			exit(0);

		if (nchar > 1 && nchar < 80)
		{

			fread(string, nchar, 1, *fp);
			string[nchar] = '\0';
			// For debugging only
			printf("\n%d\t%s", nchar, string), fflush(stdout);
			nbytes += nchar;

			if (strcmp(string, "HEADER_END") == 0)
				break;

			if (strcmp(string, "tsamp") == 0)
			{
				fread(&temp, sizeof(double), 1, *fp);
				*tsamp = (float) temp;
			}
			else if (strcmp(string, "tstart") == 0)
			{
				fread(&temp, sizeof(double), 1, *fp);
				*tstart = (float) temp;
			}
			else if (strcmp(string, "fch1") == 0)
			{
				fread(&temp, sizeof(double), 1, *fp);
				*fch1 = (float) temp;
			}
			else if (strcmp(string, "foff") == 0)
			{
				fread(&temp, sizeof(double), 1, *fp);
				*foff = (float) temp;
			}
			else if (strcmp(string, "nchans") == 0)
			{
				fread(nchans, sizeof(int), 1, *fp);
			}
			else if (strcmp(string, "nifs") == 0)
			{
				fread(nifs, sizeof(int), 1, *fp);
			}
			else if (strcmp(string, "nbits") == 0)
			{
				fread(nbits, sizeof(int), 1, *fp);
			}
			else if (strcmp(string, "nsamples") == 0)
			{
				fread(nsamples, sizeof(int), 1, *fp);
			}
		}
	}

	// Check that we are working with one IF channel
	if (*nifs != 1)
	{
		printf("\nERROR!! Can only work with one IF channel!\n");
		exit(1);
	}

	fgetpos(*fp, &file_loc);
        long data_start = ftell(*fp);
	if (fseek(*fp, 0, SEEK_END) != 0) 
 	{
		printf("\nERROR!! Failed to seek to end of data file\n");
		exit(1);
	}
        total_data = ftell(*fp);
	if (total_data == -1)
	{
		printf("\nERROR!! Failed to seek to end of data file\n");
		exit(1);
	}
	total_data -= data_start;

	if (( *nbits ) == 32)
	{
		*nsamp = (total_data/sizeof(float)) - 1;
	}
	else if (( *nbits ) == 8)
	{
		*nsamp = (total_data/sizeof(unsigned char)) - 1;
	}
	else
	{
		printf("\n\n======================= ERROR =======================\n");
		printf(" Currently this code only runs with 8 and 32 bit data\n");
		printf("\n=====================================================\n");
	}

	// Move the file pointer back to the end of the header
	fsetpos(*fp, &file_loc);

}
