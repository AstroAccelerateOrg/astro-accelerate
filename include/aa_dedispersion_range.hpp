#ifndef ASTRO_ACCELERATE_AA_DEDISPERSION_RANGE_HPP
#define ASTRO_ACCELERATE_AA_DEDISPERSION_RANGE_HPP

#include <iostream>
#include <string>
#include <stdio.h>

#include "aa_log.hpp"


namespace astroaccelerate {

class aa_dedispersion_range {
private:
	double c_dm_low;         /** Lowest value of the dispersion measure in this dedispersion range [pc*cm-3] */
	double c_dm_high;        /** Highest value of the dispersion measure in this dedispersion range [pc*cm-3] */
	double c_dm_step;        /** Step in dispersion measure [pc*cm-3] */
	double c_sampling_time;  /** Sampling time before binning [s] */
	int    c_inBin;          /** Binning factor */
	size_t c_nTimesamples;   /** Number of time-samples before binning */
	size_t c_nDMs;           /** Number of DM-trials in this dedispersion range */

public:
	/** \brief Constructor for Dedispersion_Range. */
	aa_dedispersion_range() {
		c_dm_step = 0;
		c_dm_low = 0;
		c_dm_high = 0;
		c_inBin = 0;
		c_nTimesamples = 0;
		c_nDMs = 0;
		c_sampling_time = 0;
	}

	/** \brief Constructor for Dedispersion_Range. */
	aa_dedispersion_range(
		double dm_low, 
		double dm_high, 
		double dm_step, 
		int inBin, 
		size_t nTimesamples, 
		size_t nDMs, 
		double sampling_time
	) {
		c_dm_step = dm_step;
		c_dm_low = dm_low;
		c_dm_high = dm_high;
		c_inBin = inBin;
		c_nTimesamples = nTimesamples;
		c_nDMs = nDMs;
		c_sampling_time = sampling_time;
	}

	/** \brief Method to assign or change values of an instance after construction. */
	void Assign(
		double dm_low, 
		double dm_high, 
		double dm_step, 
		int inBin, 
		size_t nTimesamples, 
		size_t nDMs, 
		double sampling_time
	) {
		c_dm_step = dm_step;
		c_dm_low = dm_low;
		c_dm_high = dm_high;
		c_inBin = inBin;
		c_nTimesamples = nTimesamples;
		c_nDMs = nDMs;
		c_sampling_time = sampling_time;
	}
	
	void print(){
		LOG(log_level::debug, 
		"  dedispersion range:" 
		+ std::to_string(c_dm_low) 
		+ "--" 
		+ std::to_string(c_dm_high) 
		+ ":" 
		+ std::to_string(c_dm_step) 
		+ "; Binning:" 
		+ std::to_string(c_inBin) 
		+ "; nTimesamples:" 
		+ std::to_string(c_nTimesamples) 
		+ "; nDMs:" 
		+ std::to_string(c_nDMs) 
		+ "; ts:" 
		+ std::to_string(c_sampling_time));
	}

	double dm_low() {return(c_dm_low);}
	double dm_high() {return(c_dm_high);}
	double dm_step() {return(c_dm_step);}
	double sampling_time() {return(c_sampling_time);}
	int    inBin() {return(c_inBin);}
	size_t nDMs() {return(c_nDMs);}
	size_t nTimesamples() {return(c_nTimesamples);}
};

} // namespace astroaccelerate
#endif // ASTRO_ACCELERATE_AA_DEDISPERSION_RANGE_HPP