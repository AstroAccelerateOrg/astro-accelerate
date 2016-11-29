#include "AstroAccelerate/UserInput.h"


namespace ska {
namespace astroaccelerate {

UserInput::UserInput()
{
	_multi_file 					= 1;
	_enable_debug 				= 0;
	_enable_analysis 			= 0;
	_enable_periodicity 	= 0;
	_enable_acceleration 	= 0;
	_output_dm 						= 0;
	_enable_zero_dm 			= 0;
	_nboots 							= -1;
	_ntrial_bins					= 0;
	_navdms								= 1;
	_narrow								= 0.001f;
	_agression						= 2.5;
	_nsearch							= 3;
	_inBin								= NULL;
	_outBin								= NULL;
	_power								= 2.0f;
	_sigma_cutoff					= 6.0f;
	_range								= 0;
	_user_dm_low					= NULL;
	_user_dm_high					= NULL;
	_user_dm_step					= NULL;
}

UserInput::~UserInput()
{
}

int 		UserInput::get_multi_file() const
{
	return _multi_file;
}

int 		UserInput::get_enable_debug() const
{
	return _enable_debug;
}

int 		UserInput::get_enable_analysis() const
{
	return _enable_analysis;
}

int 		UserInput::get_enable_periodicity() const
{
	return _enable_periodicity;
}

int 		UserInput::get_enable_acceleration() const
{
	return _enable_acceleration;
}

int 		UserInput::get_output_dm() const
{
	return _output_dm;
}

int 		UserInput::get_enable_zero_dm() const
{
	return _enable_zero_dm;
}

int 		UserInput::get_nboots() const
{
	return _nboots;
}

int 		UserInput::get_ntrial_bons() const
{
	return _ntrial_bins;
}

int 		UserInput::get_navdms() const
{
	return _navdms;
}

int 		UserInput::get_narrow() const
{
	return _narrow;
}

int 		UserInput::get_agression() const
{
	return _agression;
}

int			UserInput::get_nsearch() const
{
	return _nsearch;
}

int* 		UserInput::get_inBin() const
{
	return _inBin;
}

int* 		UserInput::get_outBin() const
{
	return _outBin;
}

float 	UserInput::get_power() const
{
	return _power;
}

float 	UserInput::get_sigma_cutoff() const
{
	return _sigma_cutoff;
}

int 		UserInput::get_range() const
{
	return _range;
}

float*	UserInput::get_user_dm_low() const
{
	return _user_dm_low;
}

float*	UserInput::get_user_dm_high() const
{
	return _user_dm_high;
}

float* 	UserInput::get_user_dm_step() const
{
	return _user_dm_step;
}

void 	UserInput::get_user_input(FILE** fp, int argc, char *argv[])
{
	/*
	FILE *fp_in = NULL;

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
		( *range ) = 0;
		while (!feof(fp_in))
		{
			//if (fscanf(fp_in, "%s", string) != 1)
			//	fprintf(stderr, "failed to read string\n");
			fscanf(fp_in, "%s", string);
			if (strcmp(string, "range") == 0)
				( *range )++;
		}
		rewind(fp_in);

		*user_dm_low = (float *) malloc(( *range ) * sizeof(float));
		*user_dm_high = (float *) malloc(( *range ) * sizeof(float));
		*user_dm_step = (float *) malloc(( *range ) * sizeof(float));
		*outBin = (int *) malloc(( *range ) * sizeof(int));
		*inBin = (int *) malloc(( *range ) * sizeof(int));

		for (i = 0; i < *range; i++)
		{
			if (fscanf(fp_in, "%s %f %f %f %d %d\n", string, &( *user_dm_low )[i], &( *user_dm_high )[i], &( *user_dm_step )[i], &( *inBin )[i], &( *outBin )[i]) !=6 )
				fprintf(stderr, "failed to read input\n");
		}

		rewind(fp_in);
		while (!feof(fp_in))
		{

			//if (fscanf(fp_in, "%s", string) != 1)
			//	fprintf(stderr, "failed to read string\n");
			fscanf(fp_in, "%s", string);
			if (strcmp(string, "debug") == 0)
				*enable_debug = 1;
			if (strcmp(string, "analysis") == 0)
				*enable_analysis = 1;
			if (strcmp(string, "periodicity") == 0)
				*enable_periodicity = 1;
			if (strcmp(string, "acceleration") == 0)
				*enable_acceleration = 1;
			if (strcmp(string, "output_dmt") == 0)
				*output_dmt = 1;
			if (strcmp(string, "zero_dm") == 0)
				*enable_zero_dm = 1;
			if (strcmp(string, "multi_file") == 0)
				*multi_file = 1;
			if (strcmp(string, "sigma_cutoff") == 0)
			{
				if (fscanf(fp_in, "%f", sigma_cutoff) != 1)
					fprintf(stderr, "failed to read sigma_cutoff\n");
			}
			if (strcmp(string, "narrow") == 0)
			{
				if (fscanf(fp_in, "%f", narrow) != 1)
					fprintf(stderr, "failed to read narrow\n");
			}
			if (strcmp(string, "wide") == 0)
			{
				if (fscanf(fp_in, "%f", wide) != 1)
					fprintf(stderr, "failed to read wide\n");
			}
			if (strcmp(string, "nboots") == 0)
			{
				if (fscanf(fp_in, "%d", nboots) != 1)
					fprintf(stderr, "failed to read nboots\n");
			}
			if (strcmp(string, "navdms") == 0)
			{
				if (fscanf(fp_in, "%d", navdms) != 1)
					fprintf(stderr, "failed to read navdms\n");
			}
			if (strcmp(string, "nwindows") == 0)
			{
				if (fscanf(fp_in, "%d", ntrial_bins) != 1)
					fprintf(stderr, "failed to read ntrial_bins\n");
			}
			if (strcmp(string, "nsearch") == 0)
			{
				if (fscanf(fp_in, "%d", nsearch) != 1)
					fprintf(stderr, "failed to read nsearch\n");
			}
			if (strcmp(string, "aggression") == 0)
			{
				if (fscanf(fp_in, "%f", aggression) != 1)
					fprintf(stderr, "failed to read aggression\n");
			}
			if (strcmp(string, "power") == 0)
			{
				if (fscanf(fp_in, "%f", power) != 1)
					fprintf(stderr, "failed to read power\n");
			}
			if (strcmp(string, "file") == 0)
			{
				//if (fscanf(fp_in, "%s", string) != 1)
				//	fprintf(stderr, "failed to read string\n");
				fscanf(fp_in, "%s", string);
				if (( *fp = fopen(string, "rb") ) == NULL)
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
	*/
}

} // namespace astroaccelerate
} // namespace ska
