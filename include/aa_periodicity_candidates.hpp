#ifndef ASTRO_ACCELERATE_AA_PERIODICITY_CANDIDATES_HPP
#define ASTRO_ACCELERATE_AA_PERIODICITY_CANDIDATES_HPP
#include <stdio.h>
#include <string.h>


namespace astroaccelerate {

/**
* \class aa_periodicity_candidates_data
* \brief Class for AstroAccelerate to internally manage processing periodicity candidates.
* \brief The user should not interact with this class for periodicity. 
* \author -
* \date -
*/
class aa_periodicity_candidates_data {
private:
	unsigned int *DM;  //index of the DM trial
	unsigned int *t;   //index of the time sample in the DM trial
	unsigned int *h;   //number of harmonics summed 
	float        *SNR; //SNR value of the candidate
	// temporary
	float        *all;
	int          c_nElements;
	
	size_t       c_nCandidates;
	int          c_rangeid;
	
	void Allocate(){
		all = (float*) malloc(c_nCandidates*c_nElements*sizeof(float));
	}
	
	bool Copy_from_GPU(float *d_all){
		if(all==NULL) return(false);
		cudaError_t err = cudaMemcpy(all, d_all, c_nCandidates*c_nElements*sizeof(float), cudaMemcpyDeviceToHost);
		if(err != cudaSuccess) return(false);
		else return(true);
	}
	
public:
	aa_periodicity_candidates_data() {
		DM  = NULL;
		t   = NULL;
		h   = NULL;
		SNR = NULL;
		all = NULL;
		c_nElements = 4;
		c_nCandidates = 0;
		c_rangeid     = 0;
	}
	
	//bool Add_candidates(unsigned int *d_DM, unsigned int *d_t, unsigned int *d_h, float *d_SNR, size_t nCandidates, int rangeid){
	//	return(false);
	//}
	
	bool Add_candidates(float *d_all, size_t nCandidates, int rangeid){
		c_nCandidates = nCandidates;
		c_rangeid     = rangeid;
		c_nElements   = 4;
		Allocate();
		bool error = Copy_from_GPU(d_all);
		return(error);
	}
	
	size_t nCandidates(){
		return(c_nCandidates);
	}
	
	bool Copy_candidates_data(float *destination){
		memcpy(destination, all, c_nCandidates*c_nElements*sizeof(float));
		return(true);
	}
	
	/** \brief Processes the periodicity inBin data group. */
	void Process(float *MSD, double dm_low, double dm_step, double sampling_time, size_t nTimesamples, float mod, int enable_greedy_postprocessing) {
		#pragma omp parallel for
		for(size_t c = 0; c < c_nCandidates; c++) {
			int harmonics = (int) all[c*c_nElements+3];
			float freq_correction = 0;
			
			if(enable_greedy_postprocessing){
				int datashift = (int) (all[c*c_nElements+3]/100.0);
				harmonics = harmonics - datashift*100;
				if(harmonics>0) {
					float bin_width = (1.0/sampling_time)/(nTimesamples*mod);
					float fraction = ((double) datashift)/((double) harmonics);
					freq_correction = bin_width*fraction;
				}
			}
			all[c*c_nElements+0] = all[c*c_nElements+0]*dm_step + dm_low; 
			all[c*c_nElements+1] = all[c*c_nElements+1]*(1.0/(sampling_time*nTimesamples*mod)) + freq_correction;
			all[c*c_nElements+2] = (all[c*c_nElements+2] - MSD[2*harmonics])/(MSD[2*harmonics+1]);
			all[c*c_nElements+3] = harmonics; 
		}
	}
	
	void Export_candidates_data_to_file(FILE *fp_out){
		fwrite(all, c_nElements*sizeof(float), c_nCandidates, fp_out);
	}
	
	~aa_periodicity_candidates_data(){
		if(all != NULL) {
			free(all);
			all = NULL;
		}
	}
};


typedef aa_periodicity_candidates_data* PSR_candidate_data_pointer;

/**
 * \class aa_periodicity_candidates aa_periodicity_candidates.hpp "include/aa_periodicity_candidates.hpp"
 * \brief Class that accumulates candidates from periodicity search.
 * \author AstroAccelerate team.
 */
class aa_periodicity_candidates {
private:
	const size_t c_nElements = 4;
	std::vector<PSR_candidate_data_pointer> candidate_data;
	
public:
	aa_periodicity_candidates(){
		
	}
	
	bool Add_Candidates(float *d_all, size_t nCandidates, int rangeid, float *MSD, double dm_low, double dm_step, double sampling_time, size_t nTimesamples, float mod, int enable_greedy_postprocessing){
		int last = 0;
		bool passed = false;
		candidate_data.push_back( (new aa_periodicity_candidates_data()) );
		last = ((int) candidate_data.size()) - 1;
		if(last >= 0){
			passed = candidate_data[last]->Add_candidates(d_all, nCandidates, rangeid);
			if(passed==true) {
				candidate_data[last]->Process(MSD, dm_low, dm_step, sampling_time, nTimesamples, mod, enable_greedy_postprocessing);
				return(true);
			}
			else {
				candidate_data[last]->~aa_periodicity_candidates_data();
				return(false);
			}
		}
		return(false);
	}
	
	size_t nCandidates(){
		int nLists = (int) candidate_data.size();
		size_t nCandidates = 0;
		for(int f = 0; f < nLists; f++){
			nCandidates = nCandidates + candidate_data[f]->nCandidates();
		}
		return(nCandidates);
	}
	
	bool Copy_candidates(float *list){
		int nLists = (int) candidate_data.size();
		std::vector<size_t> scan_candidates;
		size_t temp_sum = 0;
		scan_candidates.push_back(temp_sum);
		for(int f = 0; f < nLists; f++){
			temp_sum = temp_sum + candidate_data[f]->nCandidates();
			scan_candidates.push_back(temp_sum);
		}
		
		#pragma omp parallel for schedule(static, 1)
		for(int f = 0; f < nLists; f++){
			size_t pos = scan_candidates[f]*c_nElements;
			candidate_data[f]->Copy_candidates_data(&list[pos]);
		}
		return(true);
	}
	
	bool Export_candidates_to_file(const char *filename){
		FILE *fp_out;
		if((fp_out = fopen(filename, "wb")) == NULL) {
			LOG(log_level::error, "Error opening output file!\n");
		}
		
		int nLists = (int) candidate_data.size();
		for(int f = 0; f < nLists; f++){
			candidate_data[f]->Export_candidates_data_to_file(fp_out);
		}
		fclose(fp_out);
		return(true);
	}
	
	void Free_candidates(){
		int nLists = (int) candidate_data.size();
		for(int f = 0; f < nLists; f++){
			candidate_data[f]->~aa_periodicity_candidates_data();
		}
		candidate_data.clear();
	}
	
	~aa_periodicity_candidates(){
		Free_candidates();
	}
};

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PERIODICITY_STRATEGY_HPP
