#include "../DedispersionStrategyFile.h"

namespace astroaccelerate{

	DedispersionStrategyFile::DedispersionStrategyFile(FILE** fp
													  ,int argc
													  ,char *argv[]
													  ,DedispersionStrategy &dedispersion_strategy)
	{
		// get user input
		get_user_input(fp, argc, argv, dedispersion_strategy);

		// get file data
		get_file_data(fp, dedispersion_strategy);
	}

	DedispersionStrategyFile::~DedispersionStrategyFile()
	{
	}

	void DedispersionStrategyFile::get_user_input(FILE** fp
											 	 ,int argc
											 	 ,char *argv[]
											 	 ,DedispersionStrategy &dedispersion_strategy)
	{

		FILE *fp_in = nullptr;

		char string[100];
		int i;

		//{{{ Read in the command line parameters and open the input file

		if (argc < 2)
		{
			fprintf(stderr, "Need input file.\n");
			exit(0);
		}
		else if (argc == 2 && strcmp(argv[1], "-help") != 0)
		{
			if (( fp_in = fopen(argv[1], "r") ) == NULL)
			{
				fprintf(stderr, "Invalid input file!\n");
				exit(0);
			}
			dedispersion_strategy._range = 0;
			while (!feof(fp_in))
			{
				if ( fscanf(fp_in, "%s", string) == 0 )
				{
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
				if (strcmp(string, "range") == 0)
					dedispersion_strategy._range++;
			}
			rewind(fp_in);

			dedispersion_strategy._user_dm_low = (float *) malloc(( dedispersion_strategy._range ) * sizeof(float));
			dedispersion_strategy._user_dm_high = (float *) malloc(( dedispersion_strategy._range ) * sizeof(float));
			dedispersion_strategy._user_dm_step = (float *) malloc(( dedispersion_strategy._range ) * sizeof(float));
			dedispersion_strategy._out_bin = (int *) malloc(( dedispersion_strategy._range ) * sizeof(int));
			dedispersion_strategy._in_bin = (int *) malloc(( dedispersion_strategy._range ) * sizeof(int));

			for (i = 0; i < dedispersion_strategy._range; i++)
			{
				if (fscanf(fp_in, "%s %f %f %f %d %d\n", string, &(dedispersion_strategy._user_dm_low[i]), &(dedispersion_strategy._user_dm_high[i]), &(dedispersion_strategy._user_dm_step[i]), &(dedispersion_strategy._in_bin[i]), &(dedispersion_strategy._out_bin[i])) !=6 )
				{
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
			}

			rewind(fp_in);
			while (!feof(fp_in))
			{
				if ( fscanf(fp_in, "%s", string) == 0 )
				{
					fprintf(stderr, "failed to read input\n");
					exit(0);
				}
				if (strcmp(string, "debug") == 0)
					dedispersion_strategy._enable_debug = 1;
				if (strcmp(string, "analysis") == 0)
					dedispersion_strategy._enable_analysis = 1;
				if (strcmp(string, "periodicity") == 0)
					dedispersion_strategy._enable_periodicity = 1;
				if (strcmp(string, "acceleration") == 0)
					dedispersion_strategy._enable_acceleration = 1;
				if (strcmp(string, "output_dmt") == 0)
					dedispersion_strategy._output_dmt = 1;
				if (strcmp(string, "zero_dm") == 0)
					dedispersion_strategy._enable_zero_dm = 1;
				if (strcmp(string, "zero_dm_with_outliers") == 0)
					dedispersion_strategy._enable_zero_dm_with_outliers = 1;
				if (strcmp(string, "rfi") == 0)
					dedispersion_strategy._enable_rfi = 1;
				if (strcmp(string, "fdas_custom_fft") == 0)
					dedispersion_strategy._enable_fdas_custom_fft = 1;
				if (strcmp(string, "fdas_inbin") == 0)
					dedispersion_strategy._enable_fdas_inbin = 1;
				if (strcmp(string, "fdas_norm") == 0)
					dedispersion_strategy._enable_fdas_norm = 1;
				if (strcmp(string, "multi_file") == 0)
					dedispersion_strategy._multi_file = 1;
				if (strcmp(string, "sigma_cutoff") == 0)
				{
					if ( fscanf(fp_in, "%f", &(dedispersion_strategy._sigma_cutoff)) == 0 )
					{
						fprintf(stderr, "failed to read input\n");
						exit(0);
					}
				}
				if (strcmp(string, "sigma_constant") == 0)
				{
					if ( fscanf(fp_in, "%f", &(dedispersion_strategy._sigma_constant)) == 0 )
					{
						fprintf(stderr, "failed to read input\n");
						exit(0);
					}
				}
				if (strcmp(string, "max_boxcar_width_in_sec") == 0)
				{
					if ( fscanf(fp_in, "%f", &(dedispersion_strategy._max_boxcar_width_in_sec)) == 0 )
					{
						fprintf(stderr, "failed to read input\n");
						exit(0);
					}
				}
				if (strcmp(string, "narrow") == 0)
				{
					if ( fscanf(fp_in, "%f", &(dedispersion_strategy._narrow)) == 0 )
					{
						fprintf(stderr, "failed to read input\n");
						exit(0);
					}
				}
				if (strcmp(string, "wide") == 0)
				{
					if ( fscanf(fp_in, "%f", &(dedispersion_strategy._wide)) == 0 )
					{
						fprintf(stderr, "failed to read input\n");
						exit(0);
					}
				}
				if (strcmp(string, "nboots") == 0)
				{
					if ( fscanf(fp_in, "%d", &(dedispersion_strategy._nboots)) == 0 )
					{
						fprintf(stderr, "failed to read input\n");
						exit(0);
					}
				}
				if (strcmp(string, "navdms") == 0)
				{
					if ( fscanf(fp_in, "%d", &(dedispersion_strategy._navdms)) == 0 )
					{
						fprintf(stderr, "failed to read input\n");
						exit(0);
					}
				}
				if (strcmp(string, "nwindows") == 0)
				{
					if ( fscanf(fp_in, "%d", &(dedispersion_strategy._ntrial_bins)) == 0 )
					{
						fprintf(stderr, "failed to read input\n");
						exit(0);
					}
				}
				if (strcmp(string, "nsearch") == 0)
				{
					if ( fscanf(fp_in, "%d", &(dedispersion_strategy._nsearch)) == 0 )
					{
						fprintf(stderr, "failed to read input\n");
						exit(0);
					}
				}
				if (strcmp(string, "aggression") == 0)
				{
					if ( fscanf(fp_in, "%f", &(dedispersion_strategy._aggression)) == 0 )
					{
						fprintf(stderr, "failed to read input\n");
						exit(0);
					}
				}
				if (strcmp(string, "power") == 0)
				{
					if ( fscanf(fp_in, "%f", &(dedispersion_strategy._power)) == 0 )
					{
						fprintf(stderr, "failed to read input\n");
						exit(0);
					}
				}
				if (strcmp(string, "file") == 0)
				{
					if ( fscanf(fp_in, "%s", string) == 0 )
					{
						fprintf(stderr, "failed to read input\n");
						exit(0);
					}
					if (( *fp = fopen(string, "rb") ) == nullptr)
					{
						fprintf(stderr, "Invalid data file!\n");
						exit(0);
					}
				}
			}
		}
		else if (argc == 2 && strcmp(argv[1], "-help") == 0)
		{
			help();
		}
		else
		{
			fprintf(stderr, "Cannot recognise input, try \"./astro-accelerate -help.\"\n");
			exit(0);
		}
	}


	void DedispersionStrategyFile::get_file_data(FILE **fp
												,DedispersionStrategy &dedispersion_strategy)
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
						dedispersion_strategy._tsamp = (float) temp;
					}
					else if (strcmp(string, "tstart") == 0)
					{
						if (fread(&temp, sizeof(double), 1, *fp) != 1)
						{
							fprintf(stderr, "\nError while reading file\n");
							exit(0);
						}
						dedispersion_strategy._tstart = (float) temp;
					}
					else if (strcmp(string, "fch1") == 0)
					{
						if (fread(&temp, sizeof(double), 1, *fp) != 1)
						{
							fprintf(stderr, "\nError while reading file\n");
							exit(0);
						}
						dedispersion_strategy._fch1 = (float) temp;
					}
					else if (strcmp(string, "foff") == 0)
					{
						if (fread(&temp, sizeof(double), 1, *fp) != 1)
						{
							fprintf(stderr, "\nError while reading file\n");
							exit(0);
						}
						dedispersion_strategy._foff = (float) temp;
					}
					else if (strcmp(string, "nchans") == 0)
					{
						if (fread(&(dedispersion_strategy._nchans), sizeof(int), 1, *fp) != 1)
						{
							fprintf(stderr, "\nError while reading file\n");
							exit(0);
						}
					}
					else if (strcmp(string, "nifs") == 0)
					{
						if (fread(&(dedispersion_strategy._nifs), sizeof(int), 1, *fp) != 1)
						{
							fprintf(stderr, "\nError while reading file\n");
							exit(0);
						}
					}
					else if (strcmp(string, "nbits") == 0)
					{
						if (fread(&(dedispersion_strategy._nbits), sizeof(int), 1, *fp) != 1)
						{
							fprintf(stderr, "\nError while reading file\n");
							exit(0);
						}
					}
					else if (strcmp(string, "nsamples") == 0)
					{
						if (fread(&(dedispersion_strategy._nsamples), sizeof(int), 1, *fp) != 1)
						{
							fprintf(stderr, "\nError while reading file\n");
							exit(0);
						}
					}
				}
			}

			// Check that we are working with one IF channel
			if (dedispersion_strategy._nifs != 1)
			{
				printf("\nERROR!! Can only work with one IF channel!\n");
				exit(1);
			}

			fgetpos(*fp, &file_loc);
		/*
			if (( _nbits ) == 32)
			{
				// Allocate a tempory buffer to store a line of frequency data
				float *temp_buffer = (float *) malloc(( *nchans ) * sizeof(float));

				// Count how many time samples we have
				total_data = 0;
				while (!feof(*fp))
				{
					fread(temp_buffer, sizeof(float), ( *nchans ), *fp);
					total_data++;
				}
				*nsamp = total_data - 1;

				free(temp_buffer);
			}
			else
		*/
			 if (( dedispersion_strategy._nbits ) == 8)
			{
				// Allocate a tempory buffer to store a line of frequency data
				unsigned char *temp_buffer = (unsigned char *) malloc(( dedispersion_strategy._nchans ) * sizeof(unsigned char));

				total_data = 0;
				while (!feof(*fp))
				{
					if (((fread(temp_buffer, sizeof(unsigned char), ( dedispersion_strategy._nchans ), *fp)) != (dedispersion_strategy._nchans)) && (total_data == 0))
					{
						fprintf(stderr, "\nError while reading file\n");
						exit(0);
					}
					total_data++;
				}
				dedispersion_strategy._nsamp = total_data - 1;
				free(temp_buffer);
			}
			else if (( dedispersion_strategy._nbits ) == 4)
			{
				// Allocate a tempory buffer to store a line of frequency data
				// each byte stores 2 frequency data
				// assumption: nchans is a multiple of 2
				if ((dedispersion_strategy._nchans % 2) != 0)
				{
					printf("\nNumber of frequency channels must be a power of 2 with 4bit data\n");
					exit(0);
				}
				int nb_bytes = dedispersion_strategy._nchans/2;
				unsigned char *temp_buffer = (unsigned char *) malloc( nb_bytes * sizeof(unsigned char));
				total_data = 0;
				while (!feof(*fp))
				{
					if (((fread(temp_buffer, sizeof(unsigned char), nb_bytes, *fp)) != nb_bytes) && (total_data == 0))
					{
						fprintf(stderr, "\nError while reading file\n");
						exit(0);
					}
					total_data++;
				}
				dedispersion_strategy._nsamp = total_data - 1;
				free(temp_buffer);
			}
			else
			{
				printf("\n\n======================= ERROR ==================\n");
				printf(" Currently this code only runs with 4 and 8 bit data\n");
				printf("\n==================================================\n");
			}

			//free(string);

			// Move the file pointer back to the end of the header
			fsetpos(*fp, &file_loc);

	}
} // namespace astroaccelerate

