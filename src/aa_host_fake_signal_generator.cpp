#include <math.h>
#include <stdio.h>
#include "aa_host_fake_signal_generator.hpp"
#include "aa_filterbank_metadata.hpp"

#include "aa_log.hpp"

constexpr int MAX_VALUE = 255;
constexpr double PI = 3.14159265359;

namespace astroaccelerate {
	
	bool aa_fake_signal_generator::produce_mask(const float sigma, const int length){
		// sigma>0 is the shape parameter, g_mean>0 is the mean parameter
		mask_data.resize(length);
		float maximum = 0.0;
		float g_mean = 1.0;
		float wid = 2*sigma*pow(3*log(30),0.5);
	        float step = wid*sigma/(length);
	        for(int i = 0; i < length; i++){
	        	float x = (i + 1)*step;
	        	mask_data[i] = sqrt(sigma/(2*PI*pow(x,3)))*exp((-sigma*pow(x - g_mean,2))/(2*pow(g_mean,2)*x));
	        	if (maximum <= mask_data[i]) {
	                	maximum = mask_data[i];
	                	maximum_pos = i;
	 		}
	        }	
	        // normalization part
	       for (int i = 0; i < length; i++){
	               mask_data[i] = mask_data[i]/maximum;
	       	}
	m_ready_mask = true;
	return m_ready_mask;
	}	
	
	bool aa_fake_signal_generator::get_fake_signal(const aa_ddtr_strategy &strategy, const aa_fake_signal_metadata &m_fake, const aa_filterbank_metadata &m_signal, bool &dump_to_disk){

		if(!(strategy.ready())){
                        LOG(log_level::error, "Strategy is not ready. Please run strategy plan first.");
                        return m_ready;
                }

                produce_mask(m_fake.get_signal_sigma(), m_fake.get_signal_width());
                if(!(m_ready_mask)){
                        LOG(log_level::error, "Mask is not ready. Please produce mask first.");
                        return m_ready;
                }
	
		const int nchans = m_signal.nchans();
		int nsamples = m_signal.nsamples();
		double tsamp = m_signal.tsamp();

		if(m_fake.get_signal_start() + strategy.maxshift() > nsamples){
			LOG(log_level::error, "Fake signal is too short.");	
			LOG(log_level::error, "Asked for a signal lenght: \t" + std::to_string(nsamples));
			LOG(log_level::error, "Signal injected at: \t\t" + std::to_string(m_fake.get_signal_start()));
			LOG(log_level::error, "Maxshift needed for the plan: \t" + std::to_string(strategy.maxshift()));
			LOG(log_level::error, "Increase signal at least to: \t" + std::to_string(m_fake.get_signal_start() + strategy.maxshift()));
			return m_ready;
		}

		std::vector<float> dmshifts;
		std::vector<int> shift_index(nchans);

		signal_output.resize(nchans*nsamples);
		dmshifts = strategy.dmshifts();

		double dm_pos = m_fake.get_dm_pos();
		int width = m_fake.get_signal_width();
		int signal_start = m_fake.get_signal_start();
		if (signal_start == 0) {
			signal_start = -maximum_pos-2*strategy.maxshift();
		}

		float temp = dm_pos/tsamp;
		for (int i = 0; i < nchans; i++){
			shift_index[i] = floor(dmshifts[i]*temp);
		}
	
		int repeats = 1;
		int period = 0;
		if (m_fake.get_signal_period() > 0){
			// compute the signal period in samples
			period = m_fake.get_signal_period()/1000.0/tsamp;
			// shift the signal to a maxshift before and after the signal to obtain full signal covered with signal, i.e. heads and tails
			repeats = (nsamples + 2*strategy.maxshift())/(period) + 1;
		}

		for(int r = 0; r < repeats; r++){
			for(int i = 0; i < width; i++){
				const int time_pos = (i - maximum_pos + signal_start + period*r)*nchans;
			        for(int j = 0; j < nchans; j++){
					const int global_index = j + nchans*shift_index[j] + time_pos;
					if ( (global_index < nsamples*nchans) && (global_index > 0)){
						signal_output[global_index] = (unsigned short)(MAX_VALUE*mask_data[i]); //255;
					}
			        }
			}
		}

	if(dump_to_disk){
		  FILE *fp_out;
		  if (( fp_out = fopen("fake_signal.dat", "wb") ) == NULL) {
		    fprintf(stderr, "Error opening output file!\n");
		    exit(0);
		  }
		
		  for (int i = 0; i < nchans; i++)
		      for (int j = 0; j < nsamples; j++)
		              fprintf(fp_out,"%d %d %i\n", i, j, signal_output[j*nchans + i]);
		  fclose(fp_out);
	}

	m_ready = true;
	return m_ready;
	} // end of get_fake_signal

}
