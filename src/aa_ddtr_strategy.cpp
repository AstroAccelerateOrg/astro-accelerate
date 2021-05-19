#include "aa_ddtr_strategy.hpp"

#include "aa_device_info.hpp"

//#define MAX_OUTPUT_SIZE 8589934592
//#define MAX_OUTPUT_SIZE 12884901888
#define MAX_OUTPUT_SIZE 17179869184

namespace astroaccelerate {
  /**
   * Trivial constructor for aa_ddtr_strategy which will result in a strategy that cannot be ready. 
   */
  aa_ddtr_strategy::aa_ddtr_strategy() : m_ready(false), m_strategy_already_calculated(false), m_selected_device(NULL), m_configured_for_analysis(false), is_setup(false), m_maxshift(0), m_num_tchunks(0), m_total_ndms(0), m_max_dm(0.0), m_maxshift_high(0), m_max_ndms(0), m_power(0.0), m_enable_msd_baseline_noise(false) {
    
  }

  /**
   * Constructor for aa_ddtr_strategy that computes the strategy upon construction, and sets the ready state of the instance of the class.
   */
  aa_ddtr_strategy::aa_ddtr_strategy(const aa_ddtr_plan &plan, const aa_filterbank_metadata &metadata, const size_t &free_memory, const bool &enable_analysis, aa_device_info *selected_device) : m_ready(false), m_strategy_already_calculated(false), m_selected_device(selected_device), m_configured_for_analysis(enable_analysis), is_setup(false), m_metadata(metadata), m_maxshift(0), m_num_tchunks(0), m_total_ndms(0), m_max_dm(0.0), m_maxshift_high(0), m_max_ndms(0), m_power(plan.power()), m_enable_msd_baseline_noise(plan.enable_msd_baseline_noise()) {
    strategy(plan, free_memory, enable_analysis);
  }

  /**
   * Implementation of the function that computes the strategy.
   */
  bool aa_ddtr_strategy::strategy(const aa_ddtr_plan &plan, const size_t &free_memory, const bool &enable_analysis) {
    /**
     * This method relies on defining points when nsamps is a multiple of
     * nchans - bin on the diagonal or a fraction of it.
     */
    if(m_strategy_already_calculated) {
      LOG(log_level::notice, "Strategy already calculated.");
      return true;
    }

    LOG(log_level::notice, "Calculating strategy.");
	    
    //Part of the filterbank metadata
    const int nchans  = m_metadata.nchans();
    const int nsamp   = m_metadata.nsamples();
    const float fch1  = m_metadata.fch1();
    const float foff  = m_metadata.foff();
    const float tsamp = m_metadata.tsamp();
    const int nbits = m_metadata.nbits();
    
    // Custom bandpass normalization for zerodm filtering
    if(plan.bandpass_normalization_size() == (size_t) nchans){
        m_bandpass_normalization.resize(nchans);
        std::copy( plan.bandpass_data_pointer(), plan.bandpass_data_pointer() + nchans, m_bandpass_normalization.begin() );
        LOG(log_level::notice, "Using user defined bandpass normalization values.");
    }
    else {
        // custom bandpass normalization is not set or has wrong size => using default normalization values
        if(plan.bandpass_normalization_size()>0){
            LOG(log_level::notice, "Custom bandpass normalization must be of size equal to number of frequency channels.");
        }
        LOG(log_level::notice, "Using default bandpass normalization values.");
        m_bandpass_normalization.resize(nchans,(((1<<nbits)-1.0)/2.0));
    }
    
    
    if(!plan.range()) {
      //No user requested dm settings have been added, this is an invalid aa_ddtr_plan.
      return false;
    }
    
    //Plan requested DM settings
    const size_t range = plan.range();
    m_ndms.resize(range);
    
    m_dmshifts.resize(nchans);
    
    //Strategy set DM settings
    str_dm.resize(range);

    const size_t gpu_memory = free_memory - 104857600; // - 1073741824;
    
    const double SPDT_fraction = 3.0/4.0; // 1.0 for MSD plane profile validation
    //Calculate maxshift, the number of dms for this bin and the highest value of dm to be calculated in this bin
    if (m_power != 2.0) {
      // Calculate time independent dm shifts
      for (int c = 0; c < nchans; c++) {
	m_dmshifts[c] = 4148.741601f * ( ( 1.0 / pow(( fch1 + ( foff * c ) ), m_power) ) - ( 1.0 / pow(fch1, m_power) ) );
      }
    }
    else {
      // Calculate time independent dm shifts
      for (int c = 0; c < nchans; c++) {
	m_dmshifts[c] = (float) ( 4148.741601f * ( ( 1.0 / pow((double) ( fch1 + ( foff * c ) ), m_power) ) - ( 1.0 / pow((double) fch1, m_power) ) ) );
      }
    }
    
    for(int i = 0; i < (int)range; i++)    {
      float n;
      modff(      ( ( (int) ( (  plan.user_dm(i).high - plan.user_dm(i).low ) / plan.user_dm(i).step ) + SDIVINDM ) / SDIVINDM )     , &n); // This calculates number of SDIVINDM blocks per DM range
      m_ndms[i] = (int) ( (int) n * SDIVINDM ); // This is number of DM trial per DM range
      if (m_max_ndms < m_ndms[i])
	m_max_ndms = m_ndms[i]; // looking for maximum number of DM trials for memory allocation
      m_total_ndms = m_total_ndms + m_ndms[i];
    }
    LOG(log_level::dev_debug, "Maximum number of dm trials in any of the range steps: " + std::to_string(m_max_ndms));
    
    str_dm[0].low = plan.user_dm(0).low;                        //
    str_dm[0].high = str_dm[0].low + ( m_ndms[0] * ( plan.user_dm(0).step ) );   // Redefines DM plan to suit GPU
    str_dm[0].step = plan.user_dm(0).step;                      //
		for (size_t i = 1; i < range; i++) {
			str_dm[i].low = str_dm[i - 1].high;
			str_dm[i].high = str_dm[i].low + m_ndms[i]*plan.user_dm(i).step;
			str_dm[i].step = plan.user_dm(i).step;

			if (plan.user_dm(i - 1).inBin > 1) {
				m_maxshift = (int) ceil(((str_dm[i - 1].low + str_dm[i - 1].step*m_ndms[i - 1])*m_dmshifts[nchans - 1]) / tsamp);
				
				float modulo = (float) (SDIVINT*2*SNUMREG);
				m_maxshift = (int) (( (float) m_maxshift / (float) plan.user_dm(i - 1).inBin + modulo ) / modulo );
				m_maxshift = m_maxshift*(SDIVINT*2*SNUMREG)*plan.user_dm(i - 1).inBin;
				if (m_maxshift > m_maxshift_high) {
					m_maxshift_high = m_maxshift;
				}
			}
		}
    
		if (plan.user_dm(range - 1).inBin > 1) {
			m_maxshift = (int) ceil( ((str_dm[range - 1].low + str_dm[range - 1].step*m_ndms[range - 1])*m_dmshifts[nchans - 1]) / tsamp );

			float modulo = (float) (SDIVINT*2*SNUMREG);
			m_maxshift = (int) (((float)m_maxshift / (float)plan.user_dm(range-1).inBin + modulo) / modulo);
			m_maxshift = m_maxshift*(SDIVINT*2*SNUMREG)*plan.user_dm(range - 1).inBin;
			if (m_maxshift > m_maxshift_high) {
				m_maxshift_high = m_maxshift;
			}
		}
    
    if (m_maxshift_high == 0)    {
      m_maxshift_high = (int) ceil(( ( str_dm[range-1].low + str_dm[range-1].step * ( m_ndms[range - 1] ) ) * m_dmshifts[nchans - 1] ) / tsamp);
    }
    m_max_dm = ceil(str_dm[range-1].high);
    
    m_maxshift = ( m_maxshift_high + ( SNUMREG * 2 * SDIVINT ) );
    LOG(log_level::dev_debug, "Range:\t"+std::to_string(range-1) + " MAXSHIFT:\t" + std::to_string(m_maxshift) + " Scrunch value:\t" + std::to_string(plan.user_dm(range-1).inBin));
    LOG(log_level::dev_debug, "Maximum dispersive delay:\t "+std::to_string(m_maxshift * tsamp)+ "(s).");

    if (m_maxshift >= nsamp)    {
      LOG(log_level::error, "Your maximum DM trial exceeds the number of samples you have. Reduce your maximum DM trial.");
      return false;
    }
    
    LOG(log_level::dev_debug, "Diagonal DM:\t" + std::to_string(( tsamp * nchans * 0.0001205 * powf(( fch1 + ( foff * ( nchans / 2 ) ) ), 3.0) ) / ( -foff * nchans )));
    if (m_maxshift >= nsamp)    {
      LOG(log_level::error, "Your maximum DM trial exceeds the number of samples you have. Reduce your maximum DM trial.");
      return false;
    }
    
    /* Four cases:
     * 1) nchans < m_max_ndms & nsamp fits in GPU RAM
     * 2) nchans > m_max_ndms & nsamp fits in GPU RAM
     * 3) nchans < m_max_ndms & nsamp does not fit in GPU RAM
     * 4) nchans > m_max_ndms & nsamp does not fit in GPU RAM
     */
    
    unsigned int max_tsamps;
    // Allocate memory to store the t_processed ranges:
    m_t_processed.resize(range);
    
    if (nchans < ( m_max_ndms )) {
      // This means that we can cornerturn into the allocated output buffer
      // without increasing the memory needed
      
      // Maximum number of samples we can fit in our GPU RAM is then given by:
      size_t SPDT_memory_requirements = (enable_analysis ? (sizeof(float)*(m_max_ndms)*SPDT_fraction) : 0 );
      max_tsamps = (unsigned int) ( (gpu_memory) / ( sizeof(unsigned short)*nchans + sizeof(float)*(m_max_ndms) + SPDT_memory_requirements )); // maximum number of timesamples we can fit into GPU memory
	size_t output_size_per_tchunk = ((size_t)max_tsamps)*sizeof(float)*((size_t)m_max_ndms);
	if (output_size_per_tchunk > MAX_OUTPUT_SIZE) max_tsamps = (unsigned int) ( (MAX_OUTPUT_SIZE) / (sizeof(float)*(m_max_ndms)));
      
      // Check that we dont have an out of range maxshift:
      if ((unsigned int)( m_maxshift ) > max_tsamps)    {
	LOG(log_level::error, "The selected GPU does not have enough memory for this number of dispersion trials. Reduce maximum dm or increase the size of dm step.");
	return false;
      }
      
      // Next check to see if nsamp fits in GPU RAM:
      if ((unsigned int)nsamp < max_tsamps)    {
	// We have case 1)
	// Allocate memory to hold the values of nsamps to be processed
	unsigned long int local_t_processed;
	if (nbits == 4){
		local_t_processed = (unsigned long int) floor(( (double) ( nsamp - (m_maxshift) ) / (double) plan.user_dm(range-1).inBin ) / (double) ( SDIVINT*4*SNUMREG )); //number of timesamples per block	
		local_t_processed = local_t_processed * ( SDIVINT*4*SNUMREG ) * plan.user_dm(range-1).inBin;
	}
	else {
		local_t_processed = (unsigned long int) floor(( (double) ( nsamp - (m_maxshift) ) / (double) plan.user_dm(range-1).inBin ) / (double) ( SDIVINT*2*SNUMREG )); //number of timesamples per block
		local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * plan.user_dm(range-1).inBin;
	}
	for (size_t i = 0; i < range; i++)    {
	  m_t_processed[i].resize(1);
	  if (nbits == 4) {
		  m_t_processed[i][0] = (int)floor(( (float)(local_t_processed)/(float)plan.user_dm(i).inBin)/(float)(SDIVINT*4*SNUMREG));
		  m_t_processed[i][0] = m_t_processed[i][0]*(SDIVINT*4*SNUMREG);
	  }
	  else {
		  m_t_processed[i][0] = (int)floor(( (float)(local_t_processed)/(float)plan.user_dm(i).inBin)/(float)(SDIVINT*2*SNUMREG));
		  m_t_processed[i][0] = m_t_processed[i][0]*(SDIVINT*2*SNUMREG);
	  }
	}
	( m_num_tchunks ) = 1;
	LOG(log_level::dev_debug, "In 1");
      }
      else {
	// We have case 3)
	// Work out how many time samples we can fit into ram
	int samp_block_size = max_tsamps - ( m_maxshift );
        
	// Work out how many blocks of time samples we need to complete the processing
	// upto nsamp-maxshift
	//int num_blocks = (int) floor(( (float) nsamp - ( *maxshift ) )) / ( (float) ( samp_block_size ) ) + 1;
        
	// Find the common integer amount of samples between all bins
	int local_t_processed;
	if (nbits == 4){
		local_t_processed = (int) floor(((float)(samp_block_size)/(float)plan.user_dm(range-1).inBin)/(float)(SDIVINT*4*SNUMREG));
		local_t_processed = local_t_processed*(SDIVINT*4*SNUMREG)*plan.user_dm(range-1).inBin;
	}
	else {
		local_t_processed = (int) floor(((float)(samp_block_size)/(float)plan.user_dm(range-1).inBin)/(float)(SDIVINT*2*SNUMREG));
		local_t_processed = local_t_processed*(SDIVINT*2*SNUMREG)*plan.user_dm(range-1).inBin;
	}
        
	int num_blocks = (int) ( ((float) (nsamp - m_maxshift)) / ((float) local_t_processed) );
        
	// Work out the remaining fraction to be processed
	int remainder =  nsamp - (num_blocks*local_t_processed ) - (m_maxshift);
	if (nbits == 4){
		remainder = (int) floor((float)remainder/(float)plan.user_dm(range-1).inBin)/(float)( SDIVINT*4*SNUMREG);
		remainder = remainder*(SDIVINT*4*SNUMREG)*plan.user_dm(range-1).inBin;
	}
	else {
		remainder = (int) floor((float)remainder/(float)plan.user_dm(range-1).inBin) / (float) ( SDIVINT*2*SNUMREG );
		remainder = remainder * ( SDIVINT*2*SNUMREG ) * plan.user_dm(range-1).inBin;
	}
	int rem_block = 0;
	if(remainder>0) rem_block = 1;
        
	for (size_t i = 0; i < range; i++)    {
		// Allocate memory to hold the values of nsamps to be processed
		m_t_processed[i].resize(num_blocks + rem_block);
		// Remember the last block holds less!
		for (int j = 0; j < num_blocks; j++) {
			if (nbits == 4){
				m_t_processed[i][j] = (int)floor(((float)(local_t_processed)/(float)plan.user_dm(i).inBin)/(float)(SDIVINT*4*SNUMREG));
				m_t_processed[i][j] = m_t_processed[i][j]*(SDIVINT*4*SNUMREG);
			}
			else {
				m_t_processed[i][j] = (int) floor(((float)(local_t_processed)/(float)plan.user_dm(i).inBin)/(float)(SDIVINT*2*SNUMREG));
				m_t_processed[i][j] = m_t_processed[i][j]*(SDIVINT*2*SNUMREG);
			}
		}
		// fractional bit
		if(rem_block==1){
			if (nbits == 4){
				m_t_processed[i][num_blocks] = (int) floor(((float)(remainder)/(float)plan.user_dm(i).inBin)/(float)(SDIVINT*4*SNUMREG));
				m_t_processed[i][num_blocks] = m_t_processed[i][num_blocks]*(SDIVINT*4*SNUMREG);
			}
			else {
				m_t_processed[i][num_blocks] = (int) floor(((float)(remainder)/(float)plan.user_dm(i).inBin)/(float)(SDIVINT*2*SNUMREG));
				m_t_processed[i][num_blocks] = m_t_processed[i][num_blocks]*(SDIVINT*2*SNUMREG);
			}
		}
	}
	
	m_num_tchunks = num_blocks + rem_block;
	LOG(log_level::dev_debug, "In 3");
	LOG(log_level::dev_debug, "num_blocks:\t" + std::to_string(num_blocks));
      }
    }
    else {
      // This means that we cannot cornerturn into the allocated output buffer
      // without increasing the memory needed. Set the output buffer to be as large as the input buffer:
      
      // Maximum number of samples we can fit in our GPU RAM is then given by:
      size_t SPDT_memory_requirements = (enable_analysis ? (sizeof(float)*(m_max_ndms)*SPDT_fraction) : 0 );
      max_tsamps = (unsigned int) ( ( gpu_memory ) / ( nchans * ( sizeof(float) + sizeof(unsigned short) )+ SPDT_memory_requirements ));
	printf("gpu_memory: %zu. Maximum number of tsamp: %d. SPDT: %zu max_ndms: %d\n", gpu_memory, max_tsamps, SPDT_memory_requirements, m_max_ndms);
	size_t output_size_per_tchunk = ((size_t)max_tsamps)*sizeof(float)*((size_t)nchans);
	printf("output_sizePer_chunk: %zu\n", output_size_per_tchunk);
	if (output_size_per_tchunk > MAX_OUTPUT_SIZE) max_tsamps = (unsigned int)((MAX_OUTPUT_SIZE)/(sizeof(float)*nchans));
	printf("After check tsamp_max: %d\n", max_tsamps);

      
      // Check that we dont have an out of range maxshift:
      if (( m_maxshift ) > (int)max_tsamps) {
	LOG(log_level::error, "The selected GPU does not have enough memory for this number of dispersion trials. Reduce maximum dm or increase the size of dm step.");
	return false;
      }
      
      // Next check to see if nsamp fits in GPU RAM:
      if ((unsigned int)nsamp < max_tsamps) {
	// We have case 2)
	// Allocate memory to hold the values of nsamps to be processed
	int local_t_processed;
	if (nbits == 4){
		local_t_processed = (int) floor(((float)(nsamp - (m_maxshift))/(float) plan.user_dm(range-1).inBin)/(float)(SDIVINT*4*SNUMREG));
		local_t_processed = local_t_processed*(SDIVINT*4*SNUMREG)*plan.user_dm(range-1).inBin;
	}
	else {
		local_t_processed = (int) floor(((float)(nsamp - (m_maxshift))/(float) plan.user_dm(range-1).inBin)/(float)(SDIVINT*2*SNUMREG));
		local_t_processed = local_t_processed*(SDIVINT*2*SNUMREG)*plan.user_dm(range-1).inBin;
	}

	for (size_t i = 0; i < range; i++) {
	  m_t_processed[i].resize(1);
	  if (nbits == 4) {
		  m_t_processed[i][0] = (int) floor(((float)(local_t_processed)/(float)plan.user_dm(i).inBin)/(float)(SDIVINT*4*SNUMREG));
		  m_t_processed[i][0] = m_t_processed[i][0]*(SDIVINT*4*SNUMREG);
	  }
	  else {
		  m_t_processed[i][0] = (int) floor(((float)(local_t_processed)/(float)plan.user_dm(i).inBin)/(float)(SDIVINT*2*SNUMREG));
		  m_t_processed[i][0] = m_t_processed[i][0]*(SDIVINT*2*SNUMREG);
	  }
	}
	( m_num_tchunks ) = 1;
	LOG(log_level::dev_debug, "In 2");
      }
      else {
	// We have case 4)
	// Work out how many time samples we can fit into ram
	int samp_block_size = max_tsamps - ( m_maxshift );
        
	// Work out how many blocks of time samples we need to complete the processing
	// upto nsamp-maxshift
	//int num_blocks = (int) floor(( (float) nsamp - (float) ( *maxshift ) ) / ( (float) samp_block_size ));
        
	// Find the common integer amount of samples between all bins
	int local_t_processed;
	if (nbits == 4) {
		local_t_processed = (int) floor(( (float) ( samp_block_size ) / (float) plan.user_dm(range-1).inBin ) / (float) ( SDIVINT*4*SNUMREG ));
		local_t_processed = local_t_processed * ( SDIVINT*4*SNUMREG ) * plan.user_dm(range-1).inBin;
	}
	else {
		local_t_processed = (int) floor(( (float) ( samp_block_size ) / (float) plan.user_dm(range-1).inBin ) / (float) ( SDIVINT*2*SNUMREG ));
		local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * plan.user_dm(range-1).inBin;
	}
        
	// samp_block_size was not used to calculate remainder instead there is local_t_processed which might be different
	int num_blocks = (int) ( ((float) (nsamp - m_maxshift)) / ((float) local_t_processed) );
        
	// Work out the remaining fraction to be processed
	int remainder = nsamp - ( num_blocks * local_t_processed ) - ( m_maxshift );
	if (nbits == 4) {
		remainder = (int) floor((float) remainder / (float) plan.user_dm(range-1).inBin) / (float) ( SDIVINT*4*SNUMREG );
		remainder = remainder * ( SDIVINT*4*SNUMREG ) * plan.user_dm(range-1).inBin;
	}
	else {
		remainder = (int) floor((float) remainder / (float) plan.user_dm(range-1).inBin) / (float) ( SDIVINT*2*SNUMREG );
		remainder = remainder * ( SDIVINT*2*SNUMREG ) * plan.user_dm(range-1).inBin;
	}
	int rem_block = 0;
	if(remainder>0) rem_block = 1;
        
	for (size_t i = 0; i < range; i++)    {
		// Allocate memory to hold the values of nsamps to be processed
		m_t_processed[i].resize(num_blocks + rem_block);
		// Remember the last block holds less!
		for (int j = 0; j < num_blocks; j++) {
			if (nbits == 4) {
				m_t_processed[i][j] = (int) floor(( (float) ( local_t_processed ) / (float) plan.user_dm(i).inBin ) / (float) ( SDIVINT*4*SNUMREG ));
				m_t_processed[i][j] = m_t_processed[i][j] * ( SDIVINT*4*SNUMREG );
			}
			else {
				m_t_processed[i][j] = (int) floor(( (float) ( local_t_processed ) / (float) plan.user_dm(i).inBin ) / (float) ( SDIVINT*2*SNUMREG ));
				m_t_processed[i][j] = m_t_processed[i][j] * ( SDIVINT*2*SNUMREG );
			}
		}
		// fractional bit
		if(rem_block==1){
			if (nbits == 4) {
				m_t_processed[i][num_blocks] = (int) floor(( (float) ( remainder ) / (float) plan.user_dm(i).inBin ) / (float) ( SDIVINT*4*SNUMREG ));
				m_t_processed[i][num_blocks] = m_t_processed[i][num_blocks] * ( SDIVINT*4*SNUMREG );
			}
			else {
				m_t_processed[i][num_blocks] = (int) floor(( (float) ( remainder ) / (float) plan.user_dm(i).inBin ) / (float) ( SDIVINT*2*SNUMREG ));
				m_t_processed[i][num_blocks] = m_t_processed[i][num_blocks] * ( SDIVINT*2*SNUMREG );
			}
		}
	}
	m_num_tchunks = num_blocks + rem_block;
	LOG(log_level::dev_debug, "In 4");
      }
    }
    printf("\nNOTICE: Maxshift memory needed:\t%lu MB.", nchans * ( m_maxshift ) * sizeof(unsigned short) / 1024 / 1024);
    if (nchans < ( m_max_ndms ))    {
      printf("\nNOTICE: Output memory needed:\t%lu MB.", ( m_max_ndms ) * ( m_maxshift ) * sizeof(float) / 1024 / 1024);
    }
    else {
      printf("\nNOTICE: Output memory needed:\t%lu MB.\n", nchans * ( m_maxshift ) * sizeof(float) / 1024 / 1024);
    }

    // The memory that will be allocated on the GPU in function allocate_memory_gpu is given by gpu_inputsize + gpu_outputsize
    // gpu_inputsize
    //aa_device_info& device_info = aa_device_info::instance();
    int time_samps = m_t_processed[0][0] + m_maxshift;
    m_selected_device->request_memory((size_t) time_samps * (size_t)nchans * sizeof(unsigned short));
    // gpu_outputsize depends on nchans
    if (nchans < m_max_ndms) {
      if(!m_selected_device->request_memory((size_t)time_samps * (size_t)m_max_ndms * sizeof(float))) {
	std::cout << "ERROR:  Could not request memory." << std::endl;
	return false;
      }
    }
    else {
      if(!m_selected_device->request_memory((size_t)time_samps * (size_t)nchans * sizeof(float))) {
	std::cout << "ERROR:  Could not request memory." << std::endl;
	return false;
      }
    }
	
    
    //Strategy does not change inBin, outBin.
    //Re-assign original inBin, outBin to the strategy.
    for(size_t i = 0; i < str_dm.size(); i++) {
      str_dm.at(i).inBin = plan.user_dm(i).inBin;
      str_dm.at(i).outBin = plan.user_dm(i).outBin;
    }
    
    m_ready = true;
    m_strategy_already_calculated = true;
    return true;
  }

  /**
   * Returns the ready state of the instance of the class.
   */
  bool aa_ddtr_strategy::setup() {
    is_setup = true;
    if(is_setup && ready()) {
      return true;
    }
    else {
      return false;
    }
  }

} //namespace astroaccelerate
