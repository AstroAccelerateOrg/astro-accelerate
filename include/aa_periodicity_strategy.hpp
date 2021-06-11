#ifndef ASTRO_ACCELERATE_AA_PERIODICITY_STRATEGY_HPP
#define ASTRO_ACCELERATE_AA_PERIODICITY_STRATEGY_HPP

#include <stdio.h>
//#include <cufft.h>

#include "aa_strategy.hpp"
#include "aa_periodicity_plan.hpp"
#include "aa_dedispersion_range.hpp"
//#include "aa_periodicity_processing_chunks.hpp"
#include "aa_log.hpp"
#include "aa_params.hpp"
#include "aa_device_MSD_Configuration.hpp"
//#include "aa_device_MSD_plane_profile.hpp"

namespace astroaccelerate {

/**
* \class aa_periodicity_batch
* \brief Class for AstroAccelerate to manage processing periodicity data.
* \brief The user should not interact with this class for periodicity. Instead, the user should use aa_periodicity_plan and aa_periodicity_strategy.
* \author -
* \date -
*/
class aa_periodicity_batch {
public:
	size_t nDMs_per_batch;
	size_t nTimesamples;
	size_t DM_shift;
	MSD_Configuration MSD_conf;

	/** \brief Constructor for aa_periodicity_batch. */
	aa_periodicity_batch(size_t t_nDMs_per_batch, size_t t_DM_shift, size_t t_nTimesamples, int total_blocks) {
		nDMs_per_batch = t_nDMs_per_batch;
		nTimesamples   = t_nTimesamples;
		DM_shift       = t_DM_shift;
		
		// We divide nTimesamples by 2 because MSD is determined from power spectra (or because of FFT R->C)
		MSD_conf = *(new MSD_Configuration((nTimesamples>>1) + 1, nDMs_per_batch, 0, total_blocks));
	}

	/** \brief Destructor for aa_periodicity_batch. */
	~aa_periodicity_batch() {
		//delete MSD_conf;
	}
	
	/** \brief Method for printing member variable data for an instance. */
	void print(int id) {
		LOG(log_level::debug, "    Batch:" + std::to_string(id) + "; nDMs:" + std::to_string(nDMs_per_batch) + "; nTimesamples:" + std::to_string(nTimesamples) + "; DM shift:" + std::to_string(DM_shift) + ";");
	}
};

/**
* \class aa_periodicity_range
* \brief Class for AstroAccelerate to manage periodicity ranges.
* \brief The user should not use this class for interacting with periodicity. Instead, the user should use aa_periodicity_plan and aa_periodicity_strategy.
* \author -
* \date -
*/
class aa_periodicity_range {
public:
	aa_dedispersion_range range;
	int rangeid;
	std::vector<aa_periodicity_batch> batches;
	size_t total_MSD_blocks;

	void Create_Batches(size_t max_nDMs_in_range) {
		int nRepeats, nRest;
		total_MSD_blocks = 0;
		
		nRepeats = range.nDMs()/max_nDMs_in_range;
		nRest    = range.nDMs() - nRepeats*max_nDMs_in_range;
		for(int f=0; f<nRepeats; f++) {
			aa_periodicity_batch batch(max_nDMs_in_range, f*max_nDMs_in_range, range.nTimesamples(), 0);
			batches.push_back(batch);
			if((size_t) batch.MSD_conf.nBlocks_total > total_MSD_blocks) total_MSD_blocks = batch.MSD_conf.nBlocks_total;
		}
		if(nRest>0) {
			aa_periodicity_batch batch(nRest, nRepeats*max_nDMs_in_range, range.nTimesamples(), 0);
			batches.push_back(batch);
			if((size_t) batch.MSD_conf.nBlocks_total > total_MSD_blocks) total_MSD_blocks = batch.MSD_conf.nBlocks_total;
		}
	}

	aa_periodicity_range(aa_dedispersion_range t_range, size_t corrected_c_max_nTimesamples, size_t max_nDMs_in_memory, int t_rangeid) {
		double dm_low        = t_range.dm_low();
		double dm_high       = t_range.dm_high();
		double dm_step       = t_range.dm_step();
		int    inBin         = t_range.inBin();
		size_t nTimesamples  = t_range.nTimesamples(); 
		size_t nDMs          = t_range.nDMs(); 
		double sampling_time = t_range.sampling_time();
		
		nTimesamples  = corrected_c_max_nTimesamples/inBin;
		sampling_time = sampling_time*inBin;
		range.Assign(dm_low, dm_high, dm_step, inBin, nTimesamples, nDMs, sampling_time);
		
		size_t max_nDMs_in_range = max_nDMs_in_memory*inBin;
		rangeid = t_rangeid;
		Create_Batches(max_nDMs_in_range);
	}

	/** \brief Destructor for aa_periodicity_range. */
	~aa_periodicity_range() {
		batches.clear();
	}

	/** \brief Method for printing member variable data of an instance. */
	void print(void) {
		int size = (int)batches.size();
		int batchDM_0 = batches[0].nDMs_per_batch;
		
		range.print();
		
		if(size>1) {
			int batchDM_last = batches[size-1].nDMs_per_batch;
			LOG(log_level::debug, "  Periodicity search will run " + std::to_string(size) + " batches each containing " + std::to_string(batchDM_0) + " DM trials with tail of " + std::to_string(batchDM_last) + " DM trials");
		}
		else {
			LOG(log_level::debug, "  Periodicity search will run 1 batch containing " + std::to_string(batchDM_0) + " DM trials.");
		}
		
		float MSD_block_mem = (total_MSD_blocks*MSD_PARTIAL_SIZE*sizeof(float))/(1024.0*1024.0);
		LOG(log_level::debug, "  Total number of MSD blocks is " + std::to_string(total_MSD_blocks) + " which is " + std::to_string(MSD_block_mem) + "MB");
		//#ifdef GPU_PERIODICITY_SEARCH_DEBUG
			LOG(log_level::debug, "  > Batches:");
			for(int f=0; f<(int)batches.size(); f++) batches[f].print(f);
		//#endif
	}
};


typedef aa_periodicity_range* p_range_pointer;

/**
 * \class aa_periodicity_strategy aa_periodicity_strategy.hpp "include/aa_periodicity_strategy.hpp"
 * \brief Class that receives an aa_periodicity_plan object, and produces an aa_periodicity_strategy object.
 * \details A periodicity strategy is required for any pipeline running the periodicity component.
 * \author AstroAccelerate team.
 * \date 23 October 2018.
 */
class aa_periodicity_strategy : public aa_strategy {
private:
	// User configured values;
	float c_sigma_cutoff;                      /** User selected value of sigma cutoff. Any event with SNR below this value will not be selected as candidate */
	float c_sigma_outlier_rejection_threshold; /** User selected sigma for outlier rejection. Any value with sigma greater than this will be rejected. */
	int   c_nHarmonics;                        /** User selected number of harmonics performed by the periodicity search. */
	int   c_harmonic_sum_algorithm;            /** Flag for selection of the harmonic sum algorithm */
	int   c_candidate_selection_algorithm;               /** User selected flag to select reduction algorithm to select candidates. */
	bool  c_enable_outlier_rejection;          /** User selected flag to enable or disable outlier rejection when calculating mean and standard deviation. */
	bool  c_enable_interbinning;              /** Enable or disable interpolation for periodicity search */
	bool  c_enable_scalloping_loss_mitigation; /** Enable 5-point convolution which mitigates scalloping loss */
	bool  c_enable_spectrum_whitening;         /** Enable or disable spectrum whitening (removal of the red noise) for periodicity search */
	bool  c_pad_to_nearest_higher_pow2;        /** Whether the periodicity will pad data to the nearest power of two for the Fourier transform. True by default */
	
	size_t c_max_total_MSD_blocks;
	size_t c_max_nTimesamples;
	size_t corrected_c_max_nTimesamples;
	size_t max_nDMs;
	size_t max_nDMs_in_memory;
	size_t c_input_plane_size;
	size_t c_cuFFT_workarea_size;
	double sampling_time;
	//std::vector<aa_periodicity_range> P_ranges;
	std::vector<p_range_pointer> P_ranges;
	std::vector<aa_dedispersion_range> DDTR_ranges;
	
	bool c_ready; /** Ready state of the instance. */
	//--------------------------------------------------------
	
	void copy_plan_parameters(aa_periodicity_plan &p_plan){
		c_sigma_cutoff = p_plan.sigma_cutoff();
		c_sigma_outlier_rejection_threshold = p_plan.sigma_outlier_rejection_threshold();
		c_nHarmonics = p_plan.nHarmonics();
		c_harmonic_sum_algorithm = p_plan.harmonic_sum_algorithm();
		c_candidate_selection_algorithm = p_plan.candidate_selection_algorithm();
		c_enable_outlier_rejection = p_plan.enable_outlier_rejection();
		c_enable_interbinning = p_plan.enable_interbinning();
		c_enable_scalloping_loss_mitigation = p_plan.enable_scalloping_loss_mitigation();
		c_enable_spectrum_whitening = p_plan.enable_spectrum_whitening();
		c_pad_to_nearest_higher_pow2 = p_plan.pad_to_nearest_higher_pow2();
		
		int nRanges = p_plan.nRanges();
		for(int f=0; f<nRanges; f++){
			DDTR_ranges.push_back(p_plan.get_range(f));
		}
		c_max_nTimesamples = 0;
		max_nDMs = 0;
		size_t t_timesamples, t_DMs;
		
		if(nRanges > 0) {
			c_max_nTimesamples = DDTR_ranges[0].nTimesamples();
			max_nDMs = DDTR_ranges[0].nDMs();
			sampling_time = DDTR_ranges[0].sampling_time();
		}
		for(int f = 1; f < nRanges; f++){
			t_timesamples = DDTR_ranges[f].nTimesamples();
			t_DMs = DDTR_ranges[f].nDMs();
			if(t_timesamples > c_max_nTimesamples) c_max_nTimesamples = t_timesamples;
			if(t_DMs > max_nDMs) max_nDMs = t_DMs;
		}
		
		int nearest = (int) floorf(log2f((float) c_max_nTimesamples));
		corrected_c_max_nTimesamples = (size_t) powf(2.0, (float) nearest);
	}
	
	int Calculate_max_nDMs_in_memory(size_t memory_available, float multiple_float, float multiple_ushort) {
		size_t t_max_nDMs_in_memory, itemp;
		
		size_t memory_per_DM = ((corrected_c_max_nTimesamples + 2)*(multiple_float*sizeof(float) + multiple_ushort*sizeof(ushort)));
		LOG(log_level::debug, "   Memory required for one DM trial is " + std::to_string((float) memory_per_DM/(1024.0*1024.0)) + "MB");
		// 1 for real input real, 2 for complex output, 2 for complex cuFFT, 1 for peaks + 1 ushort
		t_max_nDMs_in_memory = (memory_available*0.98)/((corrected_c_max_nTimesamples+2)*(multiple_float*sizeof(float) + multiple_ushort*sizeof(ushort)));
		if((max_nDMs + PHS_NTHREADS) < t_max_nDMs_in_memory) { //if we can fit all DM trials from largest range into memory then we need to find nearest higher multiple of PHS_NTHREADS
			itemp = (int)(max_nDMs/PHS_NTHREADS);
			if((max_nDMs%PHS_NTHREADS)>0) itemp++;
			t_max_nDMs_in_memory = itemp*PHS_NTHREADS;
		}
		itemp = (int)(t_max_nDMs_in_memory/PHS_NTHREADS); // if we cannot fit all DM trials from largest range into memory we find nearest lower multiple of PHS_NTHREADS
		t_max_nDMs_in_memory = itemp*PHS_NTHREADS;
		
		return(t_max_nDMs_in_memory);
	}
	
	void Create_Periodicity_Plan() {
		c_max_total_MSD_blocks = 0;
		int nRanges = DDTR_ranges.size();
		for(int f = 0; f < nRanges; f++) {
			P_ranges.push_back( (new aa_periodicity_range(DDTR_ranges[f], corrected_c_max_nTimesamples, max_nDMs_in_memory, f)) );
			int new_range_index = (int) P_ranges.size() - 1;
			size_t temp = P_ranges[new_range_index]->total_MSD_blocks;
			if(temp > c_max_total_MSD_blocks) c_max_total_MSD_blocks = temp;
			
		}
	}
	
	size_t Get_max_c_cuFFT_workarea_size(){
		// This whole unfortunate function is necessary because as of CUDA11 cufftMakePlan1d may return size of the cuFFT workarea larger (up to 2x) for some input values of nDMs (number of FFTs calculated). This value is larger then the value returned for some other larger value of nDMs thus we cannot rely on getting size of the cuFFT workarea just for largest number of FFT calculated.
		//cufftHandle plan_input;
		//cufftResult cufft_error;
		size_t t_max_c_cuFFT_workarea_size = 0;
		size_t temporary_c_cuFFT_workarea_size = 0;
		
		int nRanges = (int) P_ranges.size();
		for(int r = 0; r < nRanges; r++) {
			int nBatches = (int) P_ranges[r]->batches.size();
			for(int b = 0; b < nBatches; b++) {
				size_t nTimesamples = P_ranges[r]->batches[b].nTimesamples;
				size_t nDMs         = P_ranges[r]->batches[b].nDMs_per_batch;
				
				//cufftCreate(&plan_input);
				//cufftSetAutoAllocation(plan_input, false);
				//cufft_error = cufftMakePlan1d(plan_input, nTimesamples, CUFFT_R2C, nDMs, &temporary_c_cuFFT_workarea_size);
				//if (CUFFT_SUCCESS != cufft_error){
				//	printf("CUFFT error: %d", cufft_error);
				//}
				//cufftDestroy(plan_input);
				temporary_c_cuFFT_workarea_size = nDMs*nTimesamples*8;
				
				if(temporary_c_cuFFT_workarea_size > t_max_c_cuFFT_workarea_size) {
					t_max_c_cuFFT_workarea_size = temporary_c_cuFFT_workarea_size;
				}
			}
		}
		
		return(t_max_c_cuFFT_workarea_size);
	}
	
	size_t Calculate_nSegments_for_spectrum_whitening(int min_segment_length, int max_segment_length, size_t nSamples){
		size_t nSegments = 3;
		
		int seg_length = min_segment_length;
		int seg_sum = 1 + seg_length;
		while( (size_t) (seg_sum + seg_length) < nSamples) {
			seg_sum = seg_sum + seg_length;
			seg_length = min_segment_length*log(seg_sum);
			if(seg_length > max_segment_length) seg_length = max_segment_length;
			nSegments++;
		}
		return(nSegments);
	}
	
	bool Find_Periodicity_Plan(size_t memory_available){
		size_t memory_allocated = 0, memory_for_data = 0;
		
		float multiple_float  = 5.5;
		float multiple_ushort = 2.0; 
		LOG(log_level::debug, "> Configuring memory usage:");
		memory_for_data = memory_available;
		do {
			LOG(log_level::debug, "   Memory available: " + std::to_string(memory_for_data) + " = " + std::to_string((float) memory_for_data/(1024.0*1024.0)) + "MB");
			max_nDMs_in_memory = Calculate_max_nDMs_in_memory(memory_for_data, multiple_float, multiple_ushort);
			if(max_nDMs_in_memory==0) {
				LOG(log_level::error, "Error not enough memory for periodicity search!");
				return(false);
			}
			c_input_plane_size = (corrected_c_max_nTimesamples + 2)*max_nDMs_in_memory;
			memory_allocated = c_input_plane_size*multiple_float*sizeof(float) + multiple_ushort*c_input_plane_size*sizeof(ushort);
			LOG(log_level::debug, "   Maximum number of DM trials which fit into memory is: " +  std::to_string(max_nDMs_in_memory));
			LOG(log_level::debug, "   Maximum corrected number of time samples: " +  std::to_string(corrected_c_max_nTimesamples));
			LOG(log_level::debug, "   Input plane size: " +  std::to_string((((float) c_input_plane_size*sizeof(float))/(1024.0*1024.0))) + " MB;")
			
			Create_Periodicity_Plan();
			c_cuFFT_workarea_size = Get_max_c_cuFFT_workarea_size();
			
			//-------- Additional memory for MSD ------
			int nDecimations = ((int) floorf(log2f((float)c_nHarmonics))) + 2;
			size_t additional_data_size = c_max_total_MSD_blocks*MSD_PARTIAL_SIZE*sizeof(float) + nDecimations*2*MSD_RESULTS_SIZE*sizeof(float) + c_nHarmonics*2*sizeof(float) + 2*sizeof(int);
			memory_allocated = memory_allocated + additional_data_size;
			LOG(log_level::debug, "   Memory available for the component: " + std::to_string((float) memory_available/(1024.0*1024.0)) + "MB (" + std::to_string(memory_available) + " bytes)");
			LOG(log_level::debug, "   Memory allocated by the component: " + std::to_string((float) memory_allocated/(1024.0*1024.0)) + "MB (" + std::to_string(memory_allocated) + " bytes)");
			
			//-------- Memory for spectrum whitening ------
			int max_segment_length = 250;
			int min_segment_length = 6;
			size_t nSegments = Calculate_nSegments_for_spectrum_whitening(min_segment_length, max_segment_length, corrected_c_max_nTimesamples);
			size_t memory_requirements_for_SW = nSegments*sizeof(int) + nSegments*sizeof(float)*max_nDMs_in_memory;
			LOG(log_level::debug, "   Memory required for spectral whitening: " + std::to_string((float) memory_requirements_for_SW/(1024.0*1024.0)) + "MB (" + std::to_string(memory_requirements_for_SW) + " bytes)");
			
			if(memory_allocated>memory_available) {
				LOG(log_level::warning, "--> Not enough memory for given configuration of periodicity plan. Calculating new plan...");
				memory_for_data = memory_for_data - additional_data_size;
			}
		} while(memory_allocated > memory_available);
		
		LOG(log_level::debug, "   Workarea size: " + std::to_string((float) c_cuFFT_workarea_size/(1024.0*1024.0)) + "MB");
		return(true);
	}

public:
	aa_periodicity_strategy() {
		c_sigma_cutoff = 10.0;
		c_sigma_outlier_rejection_threshold = 3.0;
		c_nHarmonics = 32;
		c_harmonic_sum_algorithm = 1;
		c_candidate_selection_algorithm = 1;
		c_enable_outlier_rejection = true;
		c_enable_interbinning = true;
		c_enable_scalloping_loss_mitigation = true;
		c_enable_spectrum_whitening = true;
		c_pad_to_nearest_higher_pow2 = false;
		
		c_max_total_MSD_blocks = 0;
		c_max_nTimesamples = 0;
		corrected_c_max_nTimesamples = 0;
		max_nDMs = 0;
		max_nDMs_in_memory = 0;
		c_input_plane_size = 0;
		c_cuFFT_workarea_size = 0;
		sampling_time = 0;
		
		c_ready = false;
	}
	
	/** \brief Constructor for aa_periodicity_strategy that sets all member variables upon construction. */
	aa_periodicity_strategy(aa_periodicity_plan &plan, size_t available_memory) {
		c_ready = false;
		LOG(log_level::debug, "------> PERIODICITY STRATEGY CONFIGURATION STARTED <------");
		copy_plan_parameters(plan);
		
		bool plan_finished = Find_Periodicity_Plan(available_memory);
		
		if(plan_finished){
			LOG(log_level::debug, "> Periodicity strategy configured");
			LOG(log_level::debug, " ");
			print_info(*this);
		}
		
		if ((c_nHarmonics > 0) && (c_sigma_outlier_rejection_threshold > 0) && plan_finished) {
			c_ready = true;
		} 
		else {
			LOG(log_level::error, "Invalid periodicity strategy parameters. Check the aa_periodicity_plan input parameters.");
			print_info(*this);
		}
		LOG(log_level::debug, "----------------------------------------------------------");
	}
	
	/** Static member function that prints member variables for a provided aa_periodicity_strategy. */
	static bool print_info(aa_periodicity_strategy &strategy) {
		LOG(log_level::dev_debug, "> Periodicity strategy information:");
		LOG(log_level::dev_debug, "  PSR - Sigma cutoff for candidate selection:\t" + std::to_string(strategy.sigma_cutoff()));
		LOG(log_level::dev_debug, "  PSR - Enable outlier rejection:\t\t" + (strategy.candidate_selection_algorithm() ? std::string("true") : std::string("false")));
		LOG(log_level::dev_debug, "  PSR - Sigma cutoff for outlier rejection:\t" + std::to_string(strategy.sigma_outlier_rejection_threshold()));
		LOG(log_level::dev_debug, "  PSR - Number of harmonics:\t\t\t" + std::to_string(strategy.nHarmonics()));
		LOG(log_level::dev_debug, "  PSR - Harmonic sum algorithm:\t\t\t" + std::to_string(strategy.harmonic_sum_algorithm()));
		LOG(log_level::dev_debug, "  PSR - Candidate selection algorithm:\t\t" + std::to_string(strategy.candidate_selection_algorithm()));
		LOG(log_level::dev_debug, "  PSR - Enable interbinning:\t\t\t" + (strategy.enable_interbinning() ? std::string("true") : std::string("false")));
		LOG(log_level::dev_debug, "  PSR - Enable scalloping loss mitigation:\t" + (strategy.enable_scalloping_loss_mitigation() ? std::string("true") : std::string("false")));
		LOG(log_level::dev_debug, "  PSR - Enable spectrum whitening:\t\t" + (strategy.enable_spectrum_whitening() ? std::string("true") : std::string("false")));
		LOG(log_level::dev_debug, "  PSR - Pad to the nearest higher power of 2:\t" + (strategy.pad_to_nearest_higher_pow2() ? std::string("true") : std::string("false")));
		LOG(log_level::dev_debug, " ");
		strategy.print_dedispersion_plan();
		LOG(log_level::dev_debug, " ");
		strategy.print_periodicity_ranges();
		
		return true;
	}
	
	/** \brief Method for printing member variable information. */
	void print_dedispersion_plan() {
		LOG(log_level::debug, "> De-dispersion plan used:");
		LOG(log_level::debug, "   sampling time: " + std::to_string(sampling_time));
		LOG(log_level::debug, "   processed time samples by de-dispersion: " + std::to_string(c_max_nTimesamples));
		LOG(log_level::debug, "   corrected time samples for FFT: " + std::to_string(corrected_c_max_nTimesamples));
		LOG(log_level::debug, "   maximum number of DM trials: " + std::to_string(max_nDMs));
		for(int f=0; f<(int)DDTR_ranges.size(); f++) {
			DDTR_ranges[f].print();
		}
	}
	
	/** \brief Method for printing member variable information. */
	void print_periodicity_ranges() {
		int nPranges = (int) P_ranges.size();
		LOG(log_level::debug, "> Periodicity execution information:");
		LOG(log_level::debug, "  Number of periodicity ranges: " + std::to_string(nPranges));
		for(int r=0;r<nPranges;r++){
			P_ranges[r]->print();
		}
	}
	
	/** \returns The ready state of the instance of the class. */
	bool ready() const {
		return c_ready;
	}
	
	/** \brief Performs any remaining setup needed.
	 * \returns Whether the operation was successful or not, at the moment this is the ready state of the instance. */
	bool setup() {
		return ready();
	}

	/** \returns The name of the module. */
	std::string name() const {
		return "periodicity_strategy";
	}



	//----------------------------- Getters ---------------------------
	
	int nRanges() const {
		return((int)DDTR_ranges.size());
	}

	aa_periodicity_range get_periodicity_range(size_t id) {
		return( *(P_ranges.at(id)) );
	}

	float sigma_cutoff() const {
		return c_sigma_cutoff;
	}

	float sigma_outlier_rejection_threshold() const {
		return c_sigma_outlier_rejection_threshold;
	}

	/** \returns A number of harmonics used in harmonic summing algorithm. */
	int nHarmonics() const {
		return c_nHarmonics;
	}

	/** \returns type of the harmonic sum algorithm used. */
	int harmonic_sum_algorithm() const {
		return c_harmonic_sum_algorithm;
	}

	/** \returns type of the reduction algorithm for candidate selection. */
	int candidate_selection_algorithm() const {
		return c_candidate_selection_algorithm;
	}

	/** \returns A boolean indicating whether the outlier rejection will be enabled, for an instance of aa_periodicity_strategy. */
	bool enable_outlier_rejection() const {
		return c_enable_outlier_rejection;
	}

	/** \returns A boolean indicating whether the interpolation will be enabled, for an instance of aa_periodicity_strategy. */
	bool enable_interbinning() const {
		return(c_enable_interbinning);
	}

	/** \returns A boolean indicating whether 5-point convolution scalloping loss mitigation is enabled. */
	bool enable_scalloping_loss_mitigation() const {
		return(c_enable_scalloping_loss_mitigation);
	}

	/** \returns A boolean indicating whether the spectrum whitening will be enabled, for an instance of aa_periodicity_strategy. */
	bool enable_spectrum_whitening() const {
		return(c_enable_spectrum_whitening);
	}

	/** \returns A boolean indicating whether the time-series will be padded to nearest power of two. */
	bool pad_to_nearest_higher_pow2() const {
		return(c_pad_to_nearest_higher_pow2);
	}
	
	/** \returns Memory block size required for periodicity search, but not total memory. */
	size_t input_plane_size(){
		return(c_input_plane_size);
	}
	
	/** \returns Memory block size required for MSD calculation. */
	size_t max_total_MSD_blocks(){
		return(c_max_total_MSD_blocks);
	}
	
	/** \returns Memory block size required for cuFFT calculation. */
	size_t cuFFT_workarea_size(){
		return(c_cuFFT_workarea_size);
	}
	
	/** \returns Memory block size required for cuFFT calculation. */
	size_t max_nTimesamples(){
		return(c_max_nTimesamples);
	}

};

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PERIODICITY_STRATEGY_HPP
