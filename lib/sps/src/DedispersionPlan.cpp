#include "../DedispersionPlan.h"

#include "../../AstroAccelerate/params.h"
#include "../../AstroAccelerate/headers_mains.h"

namespace ska {
namespace astroaccelerate {
namespace sps {

	DedispersionPlan::DedispersionPlan()
	{
		_maxshift 	 = 0;
		_dm_low 	 = NULL;
		_dm_high 	 = NULL;
		_dm_step 	 = NULL;
		_dmshifts 	 = NULL;
		_ndms 		 = NULL;
		_max_ndms 	 = 0;
		_total_ndms	 = 0;
		_max_dm		 =	0.0f;
		_range 		 = 0;
		_t_processed = NULL;
		_nbits 		 = 0;
		_nifs 		 = 0;
		_tstart 	 = 0.0f;
		_tsamp 		 = 0.0f;
		_nsamp 		 = 0;
		_nsamples 	 = 0;
		_max_samps 	 = 0;
		_nchans 	 = 0;
		_fch1 		 = 0.0f;
		_foff 		 = 0.0f;
		_num_tchunks = 0;
		_power		 = 2.0f;
	}

	DedispersionPlan::~DedispersionPlan()
	{
		// free all the pointers
		free(_dm_low);
		free(_dm_high);
		free(_dm_step);
		free(_dmshifts);
		// Probably not so trivial (int **), should free all the _t_processed[i] first
		free(_t_processed);
	}

	// Setters
	void DedispersionPlan::set_maxshift(int maxshift)
	{
		_maxshift = maxshift;
	}

	void DedispersionPlan::set_dm_low(float* dm_low)
	{
		_dm_low = dm_low;
	}

	void DedispersionPlan::set_dm_high(float* dm_high)
	{
		_dm_high = dm_high;
	}

	void DedispersionPlan::set_dm_step(float* dm_step)
	{
		_dm_step = dm_step;
	}

	void DedispersionPlan::set_dmshifts(float* dmshifts)
	{
		_dmshifts = dmshifts;
	}

	void DedispersionPlan::set_ndms(int* ndms)
	{
		_ndms = ndms;
	}

	void DedispersionPlan::set_max_ndms(int max_ndms)
	{
		_max_ndms = max_ndms;
	}

	void DedispersionPlan::set_total_ndms(int total_ndms)
	{
		_total_ndms = total_ndms;
	}

	void DedispersionPlan::set_range(int range)
	{
		_range = range;
	}

	void DedispersionPlan::set_t_processed(int** t_processed)
	{
		_t_processed = t_processed;
	}

	void DedispersionPlan::set_nbits(int nbits)
	{
		_nbits = nbits;
	}

	void DedispersionPlan::set_nifs(int nifs)
	{
		_nifs = nifs;
	}

	void DedispersionPlan::set_tstart(float tstart)
	{
		_tstart = tstart;
	}

	void DedispersionPlan::set_tsamp(float tsamp)
	{
		_tsamp = tsamp;
	}

	void DedispersionPlan::set_nsamp(int nsamp)
	{
		_nsamp = nsamp;
	}

	void DedispersionPlan::set_nsamples(int nsamples)
	{
		_nsamples = nsamples;
	}

	void DedispersionPlan::set_max_samps(int max_samps)
	{
		_max_samps = max_samps;
	}

	void DedispersionPlan::set_nchans(int nchans)
	{
		_nchans = nchans;
	}

	void DedispersionPlan::set_fch1(float fch1)
	{
		_fch1 = fch1;
	}

	void DedispersionPlan::set_foff(float foff)
	{
		_foff = foff;
	}

	void DedispersionPlan::set_num_tchunks(unsigned int num_tchunks)
	{
		_num_tchunks = num_tchunks;
	}

	void DedispersionPlan::set_power(float power)
	{
		_power = power;
	}

	// Getters

	int DedispersionPlan::get_maxshift() const
	{
		return _maxshift;
	}

	float* DedispersionPlan::get_dm_low() const
	{
		return _dm_low;
	}

	float* DedispersionPlan::get_dm_high() const
	{
		return _dm_high;
	}

	float* DedispersionPlan::get_dm_step() const
	{
		return _dm_step;
	}

	float* DedispersionPlan::get_dmshifts() const
	{
		return _dmshifts;
	}

	int* DedispersionPlan::get_ndms() const
	{
		return _ndms;
	}

	int DedispersionPlan::get_max_ndms() const
	{
		return _max_ndms;
	}

	int DedispersionPlan::get_total_ndms() const
	{
		return _total_ndms;
	}

	float DedispersionPlan::get_max_dm() const
	{
		return _max_dm;
	}
	int DedispersionPlan::get_range() const
	{
		return _range;
	}

	int** DedispersionPlan::get_t_processed() const
	{
		return _t_processed;
	}

	int DedispersionPlan::get_nbits() const
	{
		return _nbits;
	}

	int DedispersionPlan::get_nifs() const
	{
		return _nifs;
	}

	float DedispersionPlan::get_tstart() const
	{
		return _tstart;
	}

	float DedispersionPlan::get_tsamp() const
	{
		return _tsamp;
	}

	int DedispersionPlan::get_nsamp() const
	{
		return _nsamp;
	}

	int DedispersionPlan::get_nsamples() const
	{
		return _nsamples;
	}

	int DedispersionPlan::get_max_samps() const
	{
		return _max_samps;
	}

	int DedispersionPlan::get_nchans() const
	{
		return _nchans;
	}

	float DedispersionPlan::get_fch1() const
	{
		return _fch1;
	}

	float DedispersionPlan::get_foff() const
	{
		return _foff;
	}

	unsigned int DedispersionPlan::get_num_tchunks() const
	{
		return _num_tchunks;
	}

	float DedispersionPlan::get_power() const
	{
		return _power;
	}

	void DedispersionPlan::get_file_data(FILE **fp)
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
			//if (fread(&nchar, sizeof(int), 1, *fp) != 1)
			//	fprintf(stderr, "Error while reading file\n");
			fread(&nchar, sizeof(int), 1, *fp);
			if (feof(*fp))
				exit(0);

			if (nchar > 1 && nchar < 80)
			{
				fread(string, nchar, 1, *fp);
				string[nchar] = '\0';
				nbytes += nchar;

				if (strcmp(string, "HEADER_END") == 0)
					break;

				if (strcmp(string, "tsamp") == 0)
				{
					if (fread(&temp, sizeof(double), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
					_tsamp = (float) temp;
				}
				else if (strcmp(string, "tstart") == 0)
				{
					if (fread(&temp, sizeof(double), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
					_tstart = (float) temp;
				}
				else if (strcmp(string, "fch1") == 0)
				{
					if (fread(&temp, sizeof(double), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
					_fch1 = (float) temp;
				}
				else if (strcmp(string, "foff") == 0)
				{
					if (fread(&temp, sizeof(double), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
					_foff = (float) temp;
				}
				else if (strcmp(string, "nchans") == 0)
				{
					if (fread(&_nchans, sizeof(int), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
				}
				else if (strcmp(string, "nifs") == 0)
				{
					if (fread(&_nifs, sizeof(int), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
				}
				else if (strcmp(string, "nbits") == 0)
				{
					if (fread(&_nbits, sizeof(int), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
				}
				else if (strcmp(string, "nsamples") == 0)
				{
					if (fread(&_nsamples, sizeof(int), 1, *fp) != 1)
						fprintf(stderr, "Error while reading file\n");
				}
			}
		}

		// Check that we are working with one IF channel
		if (_nifs != 1)
		{
			printf("\nERROR!! Can only work with one IF channel!\n");
			exit(1);
		}

		fgetpos(*fp, &file_loc);

		if (( _nbits ) == 32)
		{
			// Allocate a tempory buffer to store a line of frequency data
			float *temp_buffer = (float *) malloc(( _nchans ) * sizeof(float));

			// Count how many time samples we have
			total_data = 0;
			while (!feof(*fp))
			{
				fread(temp_buffer, sizeof(float), ( _nchans ), *fp);
				total_data++;
			}
			_nsamp = total_data - 1;

			free(temp_buffer);
		}
		else if (( _nbits ) == 8)
		{
			// Allocate a tempory buffer to store a line of frequency data
			unsigned char *temp_buffer = (unsigned char *) malloc(( _nchans ) * sizeof(unsigned char));

			total_data = 0;
			while (!feof(*fp))
			{
				fread(temp_buffer, sizeof(unsigned char), ( _nchans ), *fp);
				total_data++;
			}
			_nsamp = total_data - 1;

			free(temp_buffer);
		}
		else
		{
			printf("\n\n======================= ERROR =======================\n");
			printf(" Currently this code only runs with 1, 8 and 32 bit data\n");
			printf("\n=====================================================\n");
		}

		free(string);

		// Move the file pointer back to the end of the header
		fsetpos(*fp, &file_loc);

	}

	void DedispersionPlan::make_strategy(float* const user_dm_low,
	                                     float* const user_dm_high,
	                                     float* const user_dm_step,
	                                     int* 	const inBin,
	                                     size_t const gpu_memory)
	{

		// This method relies on defining points when nsamps is a multiple of
		// nchans - bin on the diagonal or a fraction of it.

		int i, j, c;
		int maxshift_high = 0;

		float n;
		float fmin = ( _fch1 + ( _foff * _nchans ) );
		float fmin_pow = powf(fmin, _power);
		float fmax_pow = powf(_fch1, _power);

		_dm_low = (float *) malloc(( _range ) * sizeof(float));
		_dm_high = (float *) malloc(( _range ) * sizeof(float));
		_dm_step = (float *) malloc(( _range ) * sizeof(float));
		_ndms = (int *) malloc(( _range ) * sizeof(int));

		_dmshifts = (float *) malloc(_nchans * sizeof(float));

		//{{{ Calculate maxshift, the number of dms for this bin and
		//the highest value of dm to be calculated in this bin

		if (_power != 2.0)
		{
			// Calculate time independent dm shifts
			for (c = 0; c < _nchans; ++c)
			{
				_dmshifts[c] = 4148.741601f * ( ( 1.0 / pow(( _fch1 + ( _foff * c ) ), _power) ) - ( 1.0 / pow(_fch1, _power) ) );
			}
		}
		else
		{
			// Calculate time independent dm shifts
			for (c = 0; c < _nchans; ++c)
			{
				_dmshifts[c] = (float) ( 4148.741601f * ( ( 1.0 / pow((double) ( _fch1 + ( _foff * c ) ), _power) ) - ( 1.0 / pow((double) _fch1, _power) ) ) );
			}
		}

		for (i = 0; i < _range; ++i)
		{
			modff(( ( (int) ( ( user_dm_high[i] - user_dm_low[i] ) / user_dm_step[i] ) + SDIVINDM ) / SDIVINDM ), &n);
			_ndms[i] = (int) ( (int) n * SDIVINDM );
			if (_max_ndms < _ndms[i])
				_max_ndms = _ndms[i];
			_total_ndms = _total_ndms + _ndms[i];
		}

		_dm_low[0]  = user_dm_low[0];
		_dm_high[0] = _ndms[0] * user_dm_step[0];
		_dm_step[0] = user_dm_step[0];
		for (i = 1; i < _range; ++i)
		{
			_dm_low[i]  = _dm_high[i - 1];
			_dm_high[i] = _dm_low[i] + _ndms[i] * user_dm_step[i];
			_dm_step[i] = user_dm_step[i];

			if (inBin[i - 1] > 1)
			{
				_maxshift = (int) ceil(( ( _dm_low[i - 1] + _dm_step[i - 1] * _ndms[i - 1] ) * _dmshifts[_nchans - 1] ) / ( _tsamp ));
				_maxshift = (int) ceil((float) ( _maxshift + ( (float) ( SDIVINT * ( SNUMREG ) ) ) ) / (float) inBin[i - 1]) / (float) ( SDIVINT * ( SNUMREG ) );
				_maxshift = ( _maxshift ) * ( SDIVINT * ( SNUMREG ) ) * inBin[i - 1];
				if (( _maxshift ) > maxshift_high)
					maxshift_high = ( _maxshift );
			}
		}

		if (inBin[_range - 1] > 1)
		{
			_maxshift = (int) ceil(( ( _dm_low[_range - 1] + _dm_step[_range - 1] * _ndms[_range - 1] ) * _dmshifts[_nchans - 1] ) / ( _tsamp ));
			_maxshift = (int) ceil((float) ( _maxshift + ( (float) ( SDIVINT * ( SNUMREG ) ) ) ) / (float) inBin[_range - 1]) / (float) ( SDIVINT * ( SNUMREG ) );
			_maxshift = _maxshift * ( SDIVINT * ( SNUMREG ) ) * inBin[_range - 1];
			if (( _maxshift ) > maxshift_high)
				maxshift_high = _maxshift;
		}

		if (maxshift_high == 0)
		{
			maxshift_high = (int) ceil(( ( _dm_low[_range - 1] + _dm_step[_range - 1] * ( _ndms[_range - 1] ) ) * _dmshifts[_nchans - 1] ) / _tsamp);
		}
		_max_dm = ceil(_dm_high[_range - 1]);
		_maxshift = ( maxshift_high + 2 * ( SNUMREG * SDIVINT ) );
		if (_maxshift >= _nsamp)
		{
			printf("\n\nERROR!! Your maximum DM trial exceeds the number of samples you have.\nReduce your maximum DM trial\n\n");
			exit(1);
		}

		if (_maxshift >= _nsamp)
		{
			printf("ERROR!! Your maximum DM trial exceeds the number of samples you have.\nReduce your maximum DM trial");
			exit(1);
		}

		/* Four cases:
		 * 1) nchans < max_ndms & nsamp fits in GPU RAM
		 * 2) nchans > max_ndms & nsamp fits in GPU RAM
		 * 3) nchans < max_ndms & nsamp does not fit in GPU RAM
		 * 4) nchans > max_ndms & nsamp does not fit in GPU RAM
		 */

		int max_tsamps;

		// Allocate memory to store the t_processed ranges:
		_t_processed = (int **) malloc(_range * sizeof(int *));

		if (_nchans < ( _max_ndms ))
		{
			// This means that we can cornerturn into the allocated output buffer
			// without increasing the memory needed

			// Maximum number of samples we can fit in our GPU RAM is then given by:
			max_tsamps = (int) ( ( gpu_memory ) / ( sizeof(unsigned short) * ( _max_ndms  + _nchans ) ) );

			// Check that we dont have an out of range maxshift:
			if ( _maxshift > max_tsamps)
			{
				printf("\nERROR!! Your GPU doens't have enough memory for this number of dispersion trials.");
				printf("\nReduce your maximum dm or increase the size of your dm step");
				exit(0);
			}

			// Next check to see if nsamp fits in GPU RAM:
			if (_nsamp < max_tsamps)
			{
				// We have case 1)
				// Allocate memory to hold the values of nsamps to be processed
				int local_t_processed = (int) floor(( (float) ( _nsamp - _maxshift ) / (float) inBin[_range - 1] ) / (float) ( SDIVINT * ( SNUMREG ) ));
				local_t_processed = local_t_processed * ( SDIVINT * ( SNUMREG ) ) * inBin[_range - 1];
				for (i = 0; i < _range; ++i)
				{
					_t_processed[i]    = (int *) malloc(sizeof(int));
					_t_processed[i][0] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT * ( SNUMREG ) ));
					_t_processed[i][0] = ( _t_processed )[i][0] * ( SDIVINT * ( SNUMREG ) );
				}
				_num_tchunks = 1;
			}
			else
			{
				// We have case 3)
				// Work out how many time samples we can fit into ram
				int samp_block_size = max_tsamps - _maxshift;

				// Work out how many blocks of time samples we need to complete the processing
				// upto nsamp-maxshift
				int num_blocks = (int) floor(( (float) _nsamp - ( _maxshift ) )) / ( (float) ( samp_block_size ) ) + 1;

				// Find the common integer amount of samples between all bins
				int local_t_processed = (int) floor(( (float) ( samp_block_size ) / (float) inBin[_range - 1] ) / (float) ( SDIVINT * ( SNUMREG ) ));
				local_t_processed = local_t_processed * ( SDIVINT * ( SNUMREG ) ) * inBin[_range - 1];

				// Work out the remaining fraction to be processed
				int remainder = ( _nsamp - ( ( num_blocks - 1 ) * local_t_processed ) - ( _maxshift ) );
				remainder = (int) floor((float) remainder / (float) inBin[_range - 1]) / (float) ( SDIVINT * ( SNUMREG ) );
				remainder = remainder * ( SDIVINT * ( SNUMREG ) ) * inBin[_range - 1];

				for (i = 0; i < _range; ++i)
				{
					// Allocate memory to hold the values of nsamps to be processed
					_t_processed[i] = (int *) malloc(num_blocks * sizeof(int));
					// Remember the last block holds less!
					for (j = 0; j < num_blocks - 1; ++j)
					{
						_t_processed[i][j] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT * ( SNUMREG ) ));
						_t_processed[i][j] = _t_processed[i][j] * ( SDIVINT * ( SNUMREG ) );
					}
					// fractional bit
					_t_processed[i][num_blocks - 1] = (int) floor(( (float) ( remainder ) / (float) inBin[i] ) / (float) ( SDIVINT * ( SNUMREG ) ));
					_t_processed[i][num_blocks - 1] = _t_processed[i][num_blocks - 1] * ( SDIVINT * ( SNUMREG ) );
				}
				_num_tchunks = num_blocks;
			}
		}
		else
		{
			// This means that we cannot cornerturn into the allocated output buffer
			// without increasing the memory needed. Set the output buffer to be as large as the input buffer:

			// Maximum number of samples we can fit in our GPU RAM is then given by:
			max_tsamps = (int) ( ( gpu_memory ) / ( _nchans * ( sizeof(float) + 2 * sizeof(unsigned short) ) ) );

			// Check that we dont have an out of range maxshift:
			if (( _maxshift ) > max_tsamps)
			{
				printf("\nERROR!! Your GPU doens't have enough memory for this number of dispersion trials.");
				printf("\nReduce your maximum dm or increase the size of your dm step");
				exit(0);
			}

			// Next check to see if nsamp fits in GPU RAM:
			if (_nsamp < max_tsamps)
			{
				// We have case 2)
				// Allocate memory to hold the values of nsamps to be processed
				int local_t_processed = (int) floor(( (float) ( _nsamp - ( _maxshift ) ) / (float) inBin[_range - 1] ) / (float) ( SDIVINT * ( SNUMREG ) ));
				local_t_processed = local_t_processed * ( SDIVINT * ( SNUMREG ) ) * inBin[_range - 1];
				for (i = 0; i < _range; i++)
				{
					_t_processed[i] = (int *) malloc(sizeof(int));
					_t_processed[i][0] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT * ( SNUMREG ) ));
					_t_processed[i][0] = ( _t_processed )[i][0] * ( SDIVINT * ( SNUMREG ) );
				}
				_num_tchunks = 1;
			}
			else
			{
				// We have case 4)
				// Work out how many time samples we can fit into ram
				int samp_block_size = max_tsamps - _maxshift;

				// Work out how many blocks of time samples we need to complete the processing
				// upto nsamp-maxshift
				int num_blocks = (int) floor(( (float) _nsamp - (float) ( _maxshift ) ) / ( (float) samp_block_size ));

				// Find the common integer amount of samples between all bins
				int local_t_processed = (int) floor(( (float) ( samp_block_size ) / (float) inBin[_range - 1] ) / (float) ( SDIVINT * ( SNUMREG ) ));
				local_t_processed = local_t_processed * ( SDIVINT * ( SNUMREG ) ) * inBin[_range - 1];

				// Work out the remaining fraction to be processed
				int remainder = _nsamp - ( num_blocks * local_t_processed ) - ( _maxshift );

				for (i = 0; i < _range; ++i)
				{
					// Allocate memory to hold the values of nsamps to be processed
					_t_processed[i] = (int *) malloc(( num_blocks + 1 ) * sizeof(int));
					// Remember the last block holds less!
					for (j = 0; j < num_blocks; ++j)
					{
						_t_processed[i][j] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT * ( SNUMREG ) ));
						_t_processed[i][j] = _t_processed[i][j] * ( SDIVINT * ( SNUMREG ) );
					}
					// fractional bit
					_t_processed[i][num_blocks] = (int) floor(( (float) ( remainder ) / (float) inBin[i] ) / (float) ( SDIVINT * ( SNUMREG ) ));
					_t_processed[i][num_blocks] = _t_processed[i][num_blocks] * ( SDIVINT * ( SNUMREG ) );
				}
				_num_tchunks  = num_blocks + 1;
			}
		}
	}
} // namespace sps
} // namespace astroaccelerate
} // namespace ska
