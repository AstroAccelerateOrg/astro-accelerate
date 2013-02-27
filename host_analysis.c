#include <stdio.h>
#include <stdlib.h>

void analysis(int nsamp, int maxshift, int ndms, size_t outputsize, float *output_buffer, float dm_low, float dm_high, float dm_step, float tsamp) {

	// Calculate the total number of values
	unsigned long int vals = (nsamp-maxshift) * ndms;

	int j, dm_count;

	// Calculate the mean
	float	mean = 0;
	float	stddev = 0;
	char	filename[200];
	double	total;

	FILE *fp_out;

	total  = 0;
	for(j = 0; j < vals; j++) {
		total += output_buffer[j];
	}
	mean = total/vals;  // Mean for entire array

	// Calculate standard deviation
	total = 0;
	for(j = 0; j < vals; j++) {
		total += (output_buffer[j] - mean)*(output_buffer[j] - mean);
	}
	stddev = sqrt(total / vals); // Stddev for entire array

	// Print mean and stddev
	printf("\nMean: %f, Stddev: %f\n\n", mean, stddev);

	// Open the output file
	sprintf(filename, "output-dm_%.2f-%.2f.dat", dm_low, dm_high);
	if ((fp_out=fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}

	// Subtract dm mean from all samples and apply threshold
        for (dm_count = 0; dm_count < ndms; dm_count++) {
		for(j = 0; j < (nsamp-maxshift); j++) {
			total = output_buffer[dm_count*(nsamp-maxshift) + j] - mean;
			if (abs(total) >= (stddev * 6) ){
				fprintf(fp_out, "%lf, %f, %f\n", j * tsamp, dm_low + dm_count*dm_step, total + mean);
			}
		}
	}

	fclose(fp_out);
}
