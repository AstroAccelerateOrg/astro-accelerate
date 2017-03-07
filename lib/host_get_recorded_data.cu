#include <stdio.h>
#include <cuda_runtime.h>
//#include <omp.h>

/* Again this funtion uses a pointer to the file pointer so that it can update the position of the file pointer :)
 * Note the brakets surrounding the input_buffer, these are needed due to operator percidence...
 * (*input_buffer_odd)[(c*(*nsamp))  + total_data] = (float)temp_buffer[c];
 */

void get_recorded_data(FILE **fp, int nsamp, int nchans, int nbits, unsigned short **input_buffer, size_t *inputsize)
{

	int c;

	unsigned long int total_data;

	//{{{ Load in the raw data from the input file and transpose
	if (nbits == 8)
	{

		// Allocate a tempory buffer to store a line of frequency data
		unsigned char *temp_buffer = (unsigned char *) malloc(nchans * sizeof(unsigned char));

		// Read in the data, transpose it and store it in the input buffer
		total_data = 0;
		while (!feof(*fp))
		{

			if (fread(temp_buffer, sizeof(unsigned char), nchans, *fp) != nchans)
				break;
			for (c = 0; c < nchans; c++)
			{
				( *input_buffer )[c + total_data * ( nchans )] = (unsigned short) temp_buffer[c];
			}
			total_data++;
		}
		free(temp_buffer);
	}
	else if (nbits == 4)
	{
		// Allocate a temporary buffer to store a line of frequency data
		// each byte stores 2 frequency data
		int nb_bytes = nchans/2;
		unsigned char *temp_buffer = (unsigned char *) malloc(nb_bytes * sizeof(unsigned char));
		// Read in the data, transpose it and store it in the input buffer
		total_data = 0;
		// 00001111
		char mask = 0x0f;
		while (!feof(*fp))
		{
			if (fread(temp_buffer, sizeof(unsigned char), nb_bytes, *fp) != nb_bytes)
				break;
			for (c = 0; c < nb_bytes; c++)
			{
				// (n >> a) & ( (1 << a) - 1) -> right shift by 'a' bits, then keep the last 'b' bits
				// Here, we right shift 4 bits and keep the last 4 bits
				( *input_buffer )[ (c*2) + total_data * ( nchans )]     = (unsigned short)( (temp_buffer[c] >> 4) & mask );
				// n & ( (1 << a ) - 1)
				// Here we keep the last 4 bits
				( *input_buffer )[ (c*2) + 1 + total_data * ( nchans )] = (unsigned short)( temp_buffer[c] & mask );
			}
			total_data++;
		}
		free(temp_buffer);
	}
	else
	{
		printf("\n\n========================= ERROR =============================\n");
		printf(" This is a SKA prototype code and only runs with 4 and 8 bit data\n");
		printf("\n===============================================================\n");
	}

	//}}}
}
