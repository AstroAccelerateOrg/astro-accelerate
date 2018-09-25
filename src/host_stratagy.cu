#include "host_strategy.hpp"

DDTR_strategy::DDTR_strategy(device_DDTR_plan& ddtr_plan const&, InputDataMeta const& std::size_t gpu_memory)
                                   : m_plan(ddtr_plan)
							       , m_ready(false)
							       , m_maxshift(0)
							       , m_tchunks(0)
							       , m_max_ndms(0)
							       , m_total_ndms(ddtr_plan.total_ndms)
							       , m_max_dm(ddtr_plan.max_dm)
							       , m_power(2.0f)
							       , m_gpu_inputsize(ddtr_plan.gpu_inputsize)
							       , m_gpu_outputsize(ddtr_plan.gpu_outputsize)
							       , m_host_inputsize(ddtr_plan.host_inputsize)
							       , m_host_outputsize(ddtr_plan.host_outputsize)
							       , m_nchans(ddtr_plan.nchans)
                                   , m_dmshifts(NULL)
                                   , t_processed(NULL)
                                   , m_ndms(NULL)
                                   , m_gpu_memory(gpu_memory)
{
}

DDTR_strategy::~DDTR_strategy() {
    if(m_dmshifts!=NULL) free(m_dmshifts);
    for(int r=0; r<m_plan.ranges(); ++r){
        if(t_processed[r]!=NULL) free(t_processed[r]); 
    }
    if(t_processed!=NULL) free(t_processed);
    if(n_ndms!=NULL) free(m_ndms);
}

const InputDataMeta& DDTR_strategy::data_meta() const
{
    return m_meta;
}

void DDTR_strategy::setup() {
  strategy();
}

bool DDTR_strategy::valid() const {
  return m_ready;
}

const device_DDTR_plan& DDTR_strategy::plan() const {
  return &m_plan;
}

int DDTR_strategy::Allocate_dmshifts(){
    if(m_dmshifts) free(m_dmshifts);
    m_dmshifts = (float *) malloc( m_nchans*sizeof(float) );
    if(m_dmshifts==NULL) return(1);
    return(0);
}

int DDTR_strategy::Allocate_t_processed_outer(){
    t_processed = (size_t **) malloc( m_plan.ranges()*sizeof(size_t *) );
    if(t_processed==NULL) return(1);
    return(0);	
}

int DDTR_strategy::Allocate_t_processed_inner(size_t t_num_tchunks){
    int error = 0;
    num_tchunks = t_num_tchunks;
    for(int r=0; r<nRanges; r++){
        t_processed[r] = (size_t *) malloc( num_tchunks*sizeof(size_t) );
        if(t_processed[r] == NULL ) error++;
    }
    return(error);
}

int DDTR_strategy::totalNumberOfTimeChunks() const
{
    return(m_plan.ranges()*m_num_tchunks);
}

int DDTR_strategy::totalNumberOfDMTrials() const
{
    return((int) m_total_ndms);		
}
	
int DDTR_strategy::strategy(int *max_ndms, int *total_ndms
                            int enable_analysis)
{
	// This method relies on defining points when nsamps is a multiple of
	// m_meta.nchans - bin on the diagonal or a fraction of it.
    
	double SPDT_fraction = 3.0/4.0; // 1.0 for MSD plane profile validation
	
	int range    = m_plan->ranges();

    if(m_ndms) free(ndms);
    if(m_ndms     = (int *) malloc( range*sizeof(int) )) return 1;
	
	int i = 0;
	int j = 0;
	int c = 0;
	int m_maxshift_high = 0;

	float n = 0;
	float fmin = ( m_meta.fch1 + ( m_meta.foff * m_meta.nchans ) );
	float fmin_pow = powf(fmin, m_power);
	float fmax_pow = powf(m_meta.fch1, m_power);

    if(Allocate_dmshifts()) return 1;

	//{{{ Calculate m_maxshift, the number of dms for this bin and
	//the highest value of dm to be calculated in this bin

	if (m_power != 2.0) {
		// Calculate time independent dm shifts
		for (c = 0; c < m_m_meta.nchans; c++) {
			( *m_dmshifts )[c] = 4148.741601f * ( ( 1.0 / pow(( m_meta.fch1 + ( m_meta.foff * c ) ), m_power) ) - ( 1.0 / pow(m_meta.fch1, m_power) ) );
		}
	}
	else {
		// Calculate time independent dm shifts
		for (c = 0; c < m_meta.nchans; c++) {
			( *m_dmshifts )[c] = (float) ( 4148.741601f * ( ( 1.0 / pow((double) ( m_meta.fch1 + ( m_meta.foff * c ) ), m_power) ) - ( 1.0 / pow((double) m_meta.fch1, m_power) ) ) );
		}
	}

	for (i = 0; i < range; i++)	{
		modff(      ( ( (int) ( ( user_dm_high[i] - user_dm_low[i] ) / user_dm_step[i] ) + SDIVINDM ) / SDIVINDM )     , &n); // This calculates number of SDIVINDM blocks per DM range
		( *ndms )[i] = (int) ( (int) n * SDIVINDM ); // This is number of DM trial per DM range
		if (*max_ndms < ( *ndms )[i])
			*max_ndms = ( *ndms )[i]; // looking for maximum number of DM trials for memory allocation
		m_total_ndms += ( *ndms )[i];
	}
	printf("\nMaximum number of dm trials in any of the range steps:\t%d", *max_ndms);

	( *dm_low )[0] = user_dm_low[0];                        // 
	( *dm_high )[0] = ( *dm_low )[0] + ( ( *ndms )[0] * ( user_dm_step[0] ) );   // Redefines DM plan to suit GPU
	( *dm_step )[0] = user_dm_step[0];                      // 
	for (i = 1; i < range; i++)	{
		( *dm_low )[i] = ( *dm_high )[i - 1];
		( *dm_high )[i] = ( *dm_low )[i] + ( *ndms )[i] * user_dm_step[i];
		( *dm_step )[i] = user_dm_step[i];

		if (inBin[i - 1] > 1) {
			*m_maxshift = (int) ceil(( ( (*dm_low)[i - 1] + (*dm_step)[i - 1]*(*ndms)[i - 1] ) * ( *dmshifts )[m_meta.nchans - 1] ) / ( m_meta.tsamp ));
			*m_maxshift = (int) ceil((float) ( *m_maxshift + ( (float) ( SDIVINT*2*SNUMREG ) ) ) / (float) inBin[i - 1]) / (float) ( SDIVINT*2*SNUMREG );
			*m_maxshift = ( *m_maxshift ) * ( SDIVINT*2*SNUMREG ) * inBin[i - 1];
			if (( *m_maxshift ) > m_maxshift_high)
				m_maxshift_high = ( *m_maxshift );
		}
	}

	if (inBin[range - 1] > 1) {
		*m_maxshift = (int) ceil(( ( ( *dm_low )[range - 1] + ( *dm_step )[range - 1] * ( *ndms )[range - 1] ) * ( *dmshifts )[m_meta.nchans - 1] ) / ( m_meta.tsamp ));
		*m_maxshift = (int) ceil((float) ( *m_maxshift + ( (float) ( SDIVINT*2*SNUMREG ) ) ) / (float) inBin[range - 1]) / (float) ( SDIVINT*2*SNUMREG );
		*m_maxshift = *m_maxshift * ( SDIVINT*2*SNUMREG ) * inBin[range - 1];
		if (( *m_maxshift ) > m_maxshift_high)
			m_maxshift_high = ( *m_maxshift );
	}

	if (m_maxshift_high == 0)	{
		m_maxshift_high = (int) ceil(( ( ( *dm_low )[range - 1] + ( *dm_step )[range - 1] * ( ( *ndms )[range - 1] ) ) * ( *dmshifts )[m_meta.nchans - 1] ) / m_meta.tsamp);
	}
	m_max_dm = ceil(( *dm_high )[range - 1]);

	*m_maxshift = ( m_maxshift_high + ( SNUMREG * 2 * SDIVINT ) );
	printf("\nRange:\t%d, MAXSHIFT:\t%d, Scrunch value:\t%d", range - 1, *m_maxshift, inBin[range - 1]);
	printf("\nMaximum dispersive delay:\t%.2f (s)", *m_maxshift * m_meta.tsamp);

	if (*m_maxshift >= nsamp)	{
		printf("\n\nERROR!! Your maximum DM trial exceeds the number of samples you have.\nReduce your maximum DM trial\n\n");
		exit(1);
	}

	printf("\nDiagonal DM:\t%f", ( m_meta.tsamp * m_meta.nchans * 0.0001205 * powf(( m_meta.fch1 + ( foff * ( m_meta.nchans / 2 ) ) ), 3.0) ) / ( -foff * m_meta.nchans ));
	if (*m_maxshift >= nsamp)	{
		printf("ERROR!! Your maximum DM trial exceeds the number of samples you have.\nReduce your maximum DM trial");
		exit(1);
	}

	/* Four cases:
	 * 1) m_meta.nchans < max_ndms & nsamp fits in GPU RAM
	 * 2) m_meta.nchans > max_ndms & nsamp fits in GPU RAM
	 * 3) m_meta.nchans < max_ndms & nsamp does not fit in GPU RAM
	 * 4) m_meta.nchans > max_ndms & nsamp does not fit in GPU RAM
	 */

	unsigned int max_m_meta.tsamps;

	// Allocate memory to store the t_processed ranges:
	( *t_processed ) = (int **) malloc(range * sizeof(int *));

	if (m_meta.nchans < ( *max_ndms )) {
		// This means that we can cornerturn into the allocated output buffer 
		// without increasing the memory needed

		// Maximum number of samples we can fit in our GPU RAM is then given by:
		//max_m_meta.tsamps = (unsigned int) ( (*m_gpu_memory) / ( sizeof(unsigned short) * ( (*max_ndms) + m_meta.nchans ) ) ); // maximum number of timesamples we can fit into GPU memory
		size_t SPDT_memory_requirements = (enable_analysis==1 ? (sizeof(float)*(*max_ndms)*SPDT_fraction) : 0 );
		max_m_meta.tsamps = (unsigned int) ( (*m_gpu_memory) / ( sizeof(unsigned short)*m_meta.nchans + sizeof(float)*(*max_ndms) + SPDT_memory_requirements )); // maximum number of timesamples we can fit into GPU memory
		
		// Check that we dont have an out of range m_maxshift:
		if (( *m_maxshift ) > max_m_meta.tsamps)	{
			printf("\nERROR!! Your GPU doens't have enough memory for this number of dispersion trials.");
			printf("\nReduce your maximum dm or increase the size of your dm step");
			exit(0);
		}

		// Next check to see if nsamp fits in GPU RAM:
		if (nsamp < max_m_meta.tsamps)	{
			// We have case 1)
			// Allocate memory to hold the values of nsamps to be processed
			unsigned long int local_t_processed = (unsigned long int) floor(( (double) ( nsamp - (*m_maxshift) ) / (double) inBin[range - 1] ) / (double) ( SDIVINT*2*SNUMREG )); //number of timesamples per block
			local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * inBin[range - 1];
			for (i = 0; i < range; i++)	{
				( *t_processed )[i] = (int *) malloc(sizeof(int)); // TODO: change to size_t
				( *t_processed )[i][0] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
				( *t_processed )[i][0] = ( *t_processed )[i][0] * ( SDIVINT*2*SNUMREG );
			}
			( *num_tchunks ) = 1;
			printf("\nIn 1\n");
		}
		else {
			// We have case 3)
			// Work out how many time samples we can fit into ram 
			int samp_block_size = max_m_meta.tsamps - ( *m_maxshift );

			// Work out how many blocks of time samples we need to complete the processing
			// upto nsamp-m_maxshift
			//int num_blocks = (int) floor(( (float) nsamp - ( *m_maxshift ) )) / ( (float) ( samp_block_size ) ) + 1;

			// Find the common integer amount of samples between all bins
			int local_t_processed = (int) floor(( (float) ( samp_block_size ) / (float) inBin[range - 1] ) / (float) ( SDIVINT*2*SNUMREG ));
			local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * inBin[range - 1];
			
			int num_blocks = (int) floor(( (float) nsamp - (float)( *m_maxshift ) )) / ( (float) ( local_t_processed ) );

			// Work out the remaining fraction to be processed
			int remainder =  nsamp -  (num_blocks*local_t_processed ) - (*m_maxshift) ;
			remainder = (int) floor((float) remainder / (float) inBin[range - 1]) / (float) ( SDIVINT*2*SNUMREG );
			remainder = remainder * ( SDIVINT*2*SNUMREG ) * inBin[range - 1];

			for (i = 0; i < range; i++)	{
				// Allocate memory to hold the values of nsamps to be processed
				( *t_processed )[i] = (int *) malloc((num_blocks + 1) * sizeof(int));
				// Remember the last block holds less!
				for (j = 0; j < num_blocks; j++) {
					( *t_processed )[i][j] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
					( *t_processed )[i][j] = ( *t_processed )[i][j] * ( SDIVINT*2*SNUMREG );
				}
				// fractional bit
				( *t_processed )[i][num_blocks] = (int) floor(( (float) ( remainder ) / (float) inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
				( *t_processed )[i][num_blocks] = ( *t_processed )[i][num_blocks] * ( SDIVINT*2*SNUMREG );
			}
			( *num_tchunks ) = num_blocks + 1;
			printf("\nIn 3\n");
			printf("\nnum_blocks:\t%d", num_blocks);
		}
	}
	else {
		// This means that we cannot cornerturn into the allocated output buffer 
		// without increasing the memory needed. Set the output buffer to be as large as the input buffer:

		// Maximum number of samples we can fit in our GPU RAM is then given by:
		//max_m_meta.tsamps = (unsigned int) ( ( *m_gpu_memory ) / ( m_meta.nchans * ( sizeof(float) + 2 * sizeof(unsigned short) ) ) );
		size_t SPDT_memory_requirements = (enable_analysis==1 ? (sizeof(float)*(*max_ndms)*SPDT_fraction) : 0 );
		max_m_meta.tsamps = (unsigned int) ( ( *m_gpu_memory ) / ( m_meta.nchans * ( sizeof(float) + sizeof(unsigned short) )+ SPDT_memory_requirements ));

		// Check that we dont have an out of range m_maxshift:
		if (( *m_maxshift ) > max_m_meta.tsamps) {
			printf("\nERROR!! Your GPU doens't have enough memory for this number of dispersion trials.");
			printf("\nReduce your maximum dm or increase the size of your dm step");
			exit(0);
		}

		// Next check to see if nsamp fits in GPU RAM:
		if (nsamp < max_m_meta.tsamps) {
			// We have case 2)
			// Allocate memory to hold the values of nsamps to be processed
			int local_t_processed = (int) floor(( (float) ( nsamp - ( *m_maxshift ) ) / (float) inBin[range - 1] ) / (float) ( SDIVINT*2*SNUMREG ));
			local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * inBin[range - 1];
			for (i = 0; i < range; i++) {
				( *t_processed )[i] = (int *) malloc(sizeof(int));
				( *t_processed )[i][0] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
				( *t_processed )[i][0] = ( *t_processed )[i][0] * ( SDIVINT*2*SNUMREG );
			}
			( *num_tchunks ) = 1;
			printf("\nIn 2\n");
		}
		else {
			// We have case 4)
			// Work out how many time samples we can fit into ram 
			int samp_block_size = max_m_meta.tsamps - ( *m_maxshift );

			// Work out how many blocks of time samples we need to complete the processing
			// upto nsamp-m_maxshift
			//int num_blocks = (int) floor(( (float) nsamp - (float) ( *m_maxshift ) ) / ( (float) samp_block_size ));

			// Find the common integer amount of samples between all bins
			int local_t_processed = (int) floor(( (float) ( samp_block_size ) / (float) inBin[range - 1] ) / (float) ( SDIVINT*2*SNUMREG ));
			local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * inBin[range - 1];
			
			// samp_block_size was not used to calculate remainder instead there is local_t_processed which might be different
			int num_blocks = (int) floor(( (float) nsamp - (float) ( *m_maxshift ) ) / ( (float) local_t_processed ));

			// Work out the remaining fraction to be processed
			int remainder = nsamp - ( num_blocks * local_t_processed ) - ( *m_maxshift );
			remainder = (int) floor((float) remainder / (float) inBin[range - 1]) / (float) ( SDIVINT*2*SNUMREG );
			remainder = remainder * ( SDIVINT*2*SNUMREG ) * inBin[range - 1];

			for (i = 0; i < range; i++)	{
				// Allocate memory to hold the values of nsamps to be processed
				( *t_processed )[i] = (int *) malloc(( num_blocks + 1 ) * sizeof(int));
				// Remember the last block holds less!
				for (j = 0; j < num_blocks; j++) {
					( *t_processed )[i][j] = (int) floor(( (float) ( local_t_processed ) / (float) inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
					( *t_processed )[i][j] = ( *t_processed )[i][j] * ( SDIVINT*2*SNUMREG );
				}
				// fractional bit
				( *t_processed )[i][num_blocks] = (int) floor(( (float) ( remainder ) / (float) inBin[i] ) / (float) ( SDIVINT*2*SNUMREG ));
				( *t_processed )[i][num_blocks] = ( *t_processed )[i][num_blocks] * ( SDIVINT*2*SNUMREG );
			}
			( *num_tchunks ) = num_blocks + 1;
			printf("\nIn 4\n");
		}
	}
	printf("\nMaxshift memory needed:\t%lu MB", m_meta.nchans * ( *m_maxshift ) * sizeof(unsigned short) / 1024 / 1024);
	if (m_meta.nchans < ( *max_ndms ))	{
		printf("\nOutput memory needed:\t%lu MB", ( *max_ndms ) * ( *m_maxshift ) * sizeof(float) / 1024 / 1024);
	}
	else {
		printf("\nOutput memory needed:\t%lu MB", m_meta.nchans * ( *m_maxshift ) * sizeof(float) / 1024 / 1024);
	}

	//If method reaches this point, the plan object is in a valid state.
	m_ready = true;
}
