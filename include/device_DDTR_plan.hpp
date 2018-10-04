#ifndef ASTRO_ACCELERATE_DEVICE_DDTR_PLAN_HPP
#define ASTRO_ACCELERATE_DEVICE_DDTR_PLAN_HPP
#include <cassert>

//TODO: rename some variables. like nchans, nsamp, tsamp to something more understandable nChannels, nSamples, sampling_time
//      
class device_DDTR_plan {
  private:  
    int m_n_ranges;
    int m_i;

  public: // TODO make private
	float *user_dm_low;
	float *user_dm_high;
	float *user_dm_step;
	int *inBin;
	
	float *dm_low;
	float *dm_high;
	float *dm_step;
	
  public:
	device_DDTR_plan()
        : m_n_ranges(0)
        , m_i(0)
		, user_dm_low(NULL)
		, user_dm_high(NULL)
		, user_dm_step(NULL)
		, inBin(NULL)
		, dm_low(NULL)
		, dm_high(NULL)
		, dm_step(NULL)
    {
	}
	
	~device_DDTR_plan(){
		if(user_dm_low!=NULL) free(user_dm_low);
		if(user_dm_high!=NULL) free(user_dm_high);
		if(user_dm_step!=NULL) free(user_dm_step);
		if(inBin!=NULL) free(inBin);
		
		if(dm_low!=NULL) free(dm_low);
		if(dm_high!=NULL) free(dm_high);
		if(dm_step!=NULL) free(dm_step);
	}

    /**
     * @brief add a range to the plan
     * @param low the lowest value in the range (inclusive)
     * @param high the highest value in the range (inclusive)
     * @param step the step size to increment the DM from low up to a max ximum of high
     * @param in_bin A mystery variable unknown except to the highest of Dddtr gurus. Many years of meditation and enligtenment you must wait dear user. TODO
     * @details must call \m allocate_ranges() before calling this method
     */
    int add_range(float low, float high, float step, float in_bin)
    {
        std::assert(m_i < m_n_ranges);
        user_dm_low[m_i]  = low;
        user_dm_high[m_i] = high;
        user_dm_step[m_i] = step;
        inBin[m_i] = in_bin;
        ++m_i;
    }

    /**
     * @brief return the number of ranges that have been added
     */
	int ranges() const
    {
        std::assert(m_i == m_n_ranges - 1); // not enough ranges have been added
        return m_n_ranges;
    }

    /**
     * @brief allocate memory for the number of ranges that will be added
     * @details must be called prior to \m add_range()
     */
	int allocate_ranges(int n)
    {
        m_i = 0;
        m_n_ranges = n;
        return allocate_ranges();
    }

  private:
	int allocate_ranges()
    {
		int error=0;

        // user ranges
        if(user_dm_low) free(user_dm_low);
		user_dm_low   = (float *) malloc( m_n_ranges*sizeof(float) );
		if(user_dm_low==NULL) error++;

        if(user_dm_high) free(user_dm_high);
		user_dm_high  = (float *) malloc( m_n_ranges*sizeof(float) );
		if(user_dm_high==NULL) error++;

        if(user_dm_step) free(user_dm_step);
		user_dm_step  = (float *) malloc( m_n_ranges*sizeof(float) );
		if(user_dm_step==NULL) error++;

        if(inBin) free(inBin);
		inBin         = (int *) malloc( m_n_ranges*sizeof(int) );
		if(inBin==NULL) error++;

        // ddtr_ranges
        if(dm_low) free(dm_low);
		dm_low   = (float *) malloc( m_n_ranges*sizeof(float) );
		if(dm_low==NULL) error++;

        if(dm_high) free(dm_high);
		dm_high  = (float *) malloc( m_n_ranges*sizeof(float) );
		if(dm_high==NULL) error++;

        if(dm_step) free(dm_step);
		dm_step  = (float *) malloc( m_n_ranges*sizeof(float) );
		if(dm_step==NULL) error++;

		return(error);
	}
};

#endif
