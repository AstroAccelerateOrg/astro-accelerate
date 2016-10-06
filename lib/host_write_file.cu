/* This function writes the transformed space (dm,t) out to file in binary format.
 */

#include <stdio.h>

<<<<<<< HEAD
void write_output(int i, int t_processed, int ndms, size_t  gpu_memory, float *output_buffer, size_t gpu_outputsize, float *dm_low, float *dm_high)
{
	FILE *fp_out;
	char filename[200];
	
	sprintf(filename, "dm_%.2f-%.2f.dat", dm_low[i], dm_high[i]);
	if ((fp_out=fopen(filename, "w")) == NULL)
	{
=======
void write_output(int i, int t_processed, int ndms, size_t  gpu_memory, float *output_buffer, size_t gpu_outputsize, float *dm_low, float *dm_high) {

	FILE	*fp_out;
	char	filename[200];
	/*
	sprintf(filename, "binary-dm_%.2f-%.2f.dat", dm_low[i], dm_high[i]);
	if ((fp_out=fopen(filename, "wb")) == NULL) {
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}

	fwrite(output_buffer, sizeof(float), gpu_outputsize, fp_out);
	*/
	sprintf(filename, "dm_%.2f-%.2f.dat", dm_low[i], dm_high[i]);
	if ((fp_out=fopen(filename, "w")) == NULL) {
>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}
	
<<<<<<< HEAD
	for(int k = 0; k < ndms; k++)
	{
		for(int l = 0; l < t_processed; l++)
		{
=======
	for(int k = 0; k < ndms; k++) {
		for(int l = 0; l < t_processed; l++) {
>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229
			fprintf(fp_out, "%f ", output_buffer[k*t_processed + l]);
		}
		fprintf(fp_out, "\n");
	}
}
