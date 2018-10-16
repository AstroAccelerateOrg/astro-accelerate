#ifndef ASTRO_ACCELERATE_FILTER_BANK_DATA_META_HPP
#define ASTRO_ACCELERATE_FILTER_BANK_DATA_META_HPP

class FilterBankDataMeta {
  public:
    // members
	size_t    nsamples;  // the number of time samples
	size_t    nifs;   // the number of polarisations in each channel
	unsigned  nbits;  // the number of bits for each value
	float     tsamp;   // the sampling time
	float     tstart;  // the absolute time of the first time sample
    
    // NOTE: These have to change
    int binning_factor;
    float dm_low;
    float dm_high;
    float dm_step;
    float start_time;
    float sampling_time;

    float binned_sampling_time;

    size_t number_dms;
    size_t timesamples;

  protected:
    // we need to ensure these are not set without going through quality control
	size_t m_nchans;  // the number of frequency channels
	float  m_fch1;    // the highest frequency
	float  m_foff;    // the channel width (must be negative)

  public:
    // member functions
	FilterBankDataMeta()
		: nsamples(0)
		, nifs(1)
		, nbits(0)
		, tsamp(0)
		, tstart(0)
        , m_nchans(0)
		, m_fch1(0)
		, m_foff(0)
    {
	}
	
	~FilterBankDataMeta(){
	}

    /**
     * @brief set the description of the channels in the filter bank data
     * @param top_frequency the highest frequency of a channel (in MHz)
     * @param channel_width the width of each channel as an offset from the neighbouring channel
     *                      - must be -ve
     *                      - in MHz
     * @param number_of_channels the number of channels
     */
    void set_channels(float top_frequency, float channel_width, size_t number_of_channels)
    {
        std::assert(channel_width < 0.0);
        std::assert(top_frequency > 0.0);
        std::assert(top_frequency + number_of_channels * channel_width > 0.0);
        m_fch1 = top_frequency;
        m_foff = channel_width;
        m_nchans = number_of_channels;
    }

    inline float fch1() const { return m_fch1; }
    inline float foff() const { return m_foff; }
    inline size_t nchans() const { return m_nchans; }

    /**
     * @brief calculate the expected size of the data
     * @details warning this may be truncated if this value is > than can be contained in a size_t
     */
    size_t data_size_in_bytes() const
    {

        return nsamples * nchans * nifs * (nbits/8);
    }

    void Setup() {
        binned_sampling_time = sampling_time * binning_factor;
    }

    void SetBinningFactor(int ibin) {
        binning_factor = ibin;
    }
};

#endif // ASTRO_ACCELERATE_FILTER_BANK_DATA_META_HPP
