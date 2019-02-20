#include <math.h>
#include <stdio.h>
#include "aa_host_fake_signal_generator.hpp"

#define MAX_VALUE 255

namespace astroaccelerate {

	void aa_fake_signal_generator::produce_mask(float sigma, int length){
		mask_data.resize(length);
		float maximum = 0.0f;
		float wid = 2*sigma*pow(3*log(30),0.5);
	        float step = wid*sigma/(length);
	        for(int i = 0; i < length; i++){
	        	float x = (i + 1)*step;
	        	mask_data[i] = pow(1.0*sigma/(2*3.141592*pow(x,3)),0.5)*exp((-sigma*pow(x - 1.0,2))/(2*pow(1.0,2)*x));
	        	if (maximum <= mask_data[i]) {
	                	maximum = mask_data[i];
	                	maximum_pos = i;
	 		}
	        }	
	        // normalization part
	       for (int i = 0; i < length; i++){
	               mask_data[i] = mask_data[i]/maximum;
	       	}
	}	

	void aa_fake_signal_generator::get_fake_signal(const aa_ddtr_strategy &strategy, const aa_fake_signal_metadata &m_fake, const aa_filterbank_metadata &m_signal){

		std::vector<float> dmshifts;
		std::vector<int> shift_index;
		int nchans = m_signal.nchans();
		int nsamples = m_signal.nsamples();
		double tsamp = m_signal.tsamp();

		shift_index.resize(nchans);
		signal_output.resize(nchans*nsamples);
		dmshifts = strategy.dmshifts();

		double dm_pos = m_fake.get_dm_pos();
		int width = m_fake.get_signal_width();
		int signal_start = m_fake.get_signal_start();

		for (int i = 0; i < nchans; i++)
			shift_index[i] = floor(dmshifts[i]*dm_pos/tsamp);

		for(int i = 0; i < width; i++){
			int time_pos = (i - maximum_pos + signal_start)*nchans;
		        for(int j = 0; j < nchans; j++){
				signal_output[j + nchans*shift_index[j] + time_pos] = (unsigned short)(MAX_VALUE*mask_data[i]); //255;
		        }
		}

	}


}
