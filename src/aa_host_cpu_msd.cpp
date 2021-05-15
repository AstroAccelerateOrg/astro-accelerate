#include "aa_host_cpu_msd.hpp"
#include "aa_log.hpp"

namespace astroaccelerate {

//void textbook_msd_parallel(unsigned int first, unsigned int second, unsigned short* in_data, float* mean, double* sd, bool outlier_rejection, double old_mean, double old_stdev, double sigma){
void textbook_msd_parallel(size_t first, size_t second, unsigned short const*const in_data, float* mean){
        double sum = 0;
        //double temp_sd = 0;
        size_t nElements = 0;

        //double low  = old_mean - sigma*old_stdev;
        //double high = old_mean + sigma*old_stdev;

        #pragma omp parallel for reduction (+: sum, nElements) // +:temp_sd
        for(size_t i = 0; i < first; i++){
                sum = 0.0;
                //temp_sd = 0.0;
                nElements = 0;
                for(size_t j = 0; j < second; j++){
//                        if (outlier_rejection && (in_data[i*second + j]<low || in_data[i*second + j]>high)){
//                        } else {
                                //temp_sd += in_data[i*second + j]*in_data[i*second + j];
                                sum += in_data[i*second + j];
                                nElements +=1;
//                        }
                }
                mean[i] = sum/nElements;
                //temp_sd = temp_sd - (sum*sum)/nElements;
		if (i == 0) printf("\n\t %f %zu\n", mean[i], nElements);
                //sd[i] = sqrt(temp_sd/nElements);
        }
}


 void call_cpu_msd(float *h_new_bandpass, unsigned short const*const input_buffer, size_t nchans, size_t nsamples){
	  LOG(log_level::debug, "Performing CPU MSD...");
	  textbook_msd_parallel(nchans, nsamples, input_buffer, h_new_bandpass);
 }
} //namespace astroaccelerate
