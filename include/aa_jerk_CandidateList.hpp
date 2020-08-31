#ifndef ASTRO_ACCELERATE_A_JERK_CANDIDATELIST_HPP
#define ASTRO_ACCELERATE_A_JERK_CANDIDATELIST_HPP

namespace astroaccelerate {
	
	struct JERK_Candidate {
		float z;
		float r; //!!
		float snr;
		float w;
	};

	class JERK_CandidateSubList;
	class JERK_CandidateList;

	class JERK_CandidateSubList { // time-DM data type
	friend class JERK_CandidateList;
	private:
		int nCandidates;
		int nHarmonics;
		
	public:
		// data
		float *candidates;
		float *MSD;
		
		// Coordinate transformation [indexes] -> [time, DM]
		float w;
		float DM;
		int nFilters_z_half;
		int nFilters_z;
		float z_max_search_limit;
		float z_search_step;
		double sampling_time;
		int nTimesaples_time_dom;
		int inBin;
		
		int getNumberOfCandidates(void) { return(nCandidates); }
		int getNumberOfHarmonics(void) { return(nHarmonics); }
		
		void SetNumberOfCandidates(int t_nCandidates){
			nCandidates = t_nCandidates;
		}
		
		void Print(){
			printf("w=%f, DM=%f; sampling time=%f; inBin=%d; nTimesaples_time_dom=%d;\n", w, DM, sampling_time, inBin, nTimesaples_time_dom);
			printf("nCandidates=%d; nHarmonics=%d;\n", nCandidates, nHarmonics);
		}
		
		int AllocateMemory(){
			if(nCandidates>0){
				candidates = (float *) malloc(nCandidates*4*sizeof(float) );
				if(candidates==NULL) return(2);
				else return(0);
			}
			return(1);
		}
		
		int exportToFile(void){
			char filename[200];
			
			sprintf(filename, "JERK_sublist-w_%.2f-dm_%.2f.dat", w, DM);
			
			FILE *fp_out;
			
			if(nCandidates>0){
				if (( fp_out = fopen(filename, "wb") ) == NULL) return(1);
				fwrite(candidates, nCandidates*sizeof(float), 4, fp_out);
				fclose(fp_out);
			}
			return(0);
		}
		
		void process_candidates(){
			//printf("nCandidates:%d; sampling time: %d; inbin: %d; nTimesamples: %d; factor: %f\n", nCandidates, sampling_time, inBin, nTimesaples_time_dom, (sampling_time*((double) inBin)*((double) nTimesaples_time_dom)) );
			for(int f=0; f<nCandidates; f++){
				//printf("Candidate [%d] = [%f; %f; %f; %f]; \n", f, candidates[4*f], candidates[4*f+1], candidates[4*f+2], candidates[4*f+3]);
				candidates[4*f]     = candidates[4*f]*z_search_step - z_max_search_limit;
				candidates[4*f + 1] = candidates[4*f + 1]/(sampling_time*((double) inBin)*((double) nTimesaples_time_dom));
				candidates[4*f + 2] = candidates[4*f + 2]; // SNR
				candidates[4*f + 3] = candidates[4*f + 3]; // w
				//printf("Candidate [%d] scaled = [%f; %f; %f; %f]; \n", f, candidates[4*f], candidates[4*f+1], candidates[4*f+2], candidates[4*f+3]);
			}
		}
		
		JERK_CandidateSubList(){
			nCandidates = 0;
			candidates = NULL;
			nHarmonics = 0;
			MSD = NULL;
			w = 0;
		}

		JERK_CandidateSubList(int t_nCandidates, float *candidates_input, float t_w, float t_DM, int t_nTimesamples_time_dom, int t_nFilters_z, float t_z_max_search_limit, float t_z_search_step, float t_sampling_time, int t_inBin){
			
			//nCandidates = t_nCandidates;
			//if(nCandidates>0){
			//	candidates = (float *) malloc(nCandidates*4*sizeof(float) );
			//	if(candidates!=NULL){
			//		memcpy(candidates, candidates_input, nCandidates*4*sizeof(float) );
			//	}
			//}
			
			//----> 
			nHarmonics = 0;
			MSD = NULL;
			
			w  = t_w;
			DM = t_DM;
			
			nFilters_z = t_nFilters_z;
			nFilters_z_half = (nFilters_z-1)/2;
			z_max_search_limit = t_z_max_search_limit;
			z_search_step = t_z_search_step;
			
			sampling_time = t_sampling_time;
			inBin = t_inBin;
			nTimesaples_time_dom = t_nTimesamples_time_dom;
			//------------------------------------------------<
			
			nCandidates = t_nCandidates;
			candidates = NULL;
			if(nCandidates>0){
				candidates = (float *) malloc(nCandidates*4*sizeof(float) );
				if(candidates!=NULL){
					if( cudaSuccess != cudaMemcpy(candidates,  candidates_input,  nCandidates*4*sizeof(float), cudaMemcpyDeviceToHost) ) {
						printf("CUDA API failure!\n");
					}
				}
				process_candidates();
			}
		}
		
		~JERK_CandidateSubList(){
			if(candidates!=NULL) free(candidates);
			if(MSD!=NULL) free(MSD);
		}
	};

	typedef JERK_CandidateSubList* JERKCandidateSubListPointer;

	class JERK_CandidateList {
	private:
		std::vector<JERKCandidateSubListPointer> clp;
		
	public:
		float *allcandidates;
		size_t nCandidates;

		size_t getNumberOfCandidates(){
			return(nCandidates);
		}
		
		size_t getNumberOfSubLists(){
			return(clp.size());
		}
		
		void AddSubListFromGPU(int t_nCandidates, float *d_candidate_input, float w, float DM, int nTimesaples_time_dom, int nFilters_z, float z_max_search_limit, float z_search_step, float sampling_time, int inBin){
			clp.push_back(new JERK_CandidateSubList(t_nCandidates, d_candidate_input, w, DM, nTimesaples_time_dom, nFilters_z, z_max_search_limit, z_search_step, sampling_time, inBin));
			int pos = ((int) clp.size()) - 1;
			//clp[pos]->Print();
		}
		
		JERK_CandidateSubList* getSubList(size_t index){
			if(index<clp.size())
				return(clp[index]);
			return(NULL);
		}
		
		void exportSubLists(){
			for(size_t f=0; f<clp.size(); f++){
				clp[f]->exportToFile();
			}
		}
		
		int poolCandidates(){
			size_t nSubLists = clp.size();
			nCandidates = 0;
			for(size_t f=0; f<nSubLists; f++){
				nCandidates = nCandidates + clp[f]->getNumberOfCandidates();
			}
			allcandidates = (float *) malloc(nCandidates*4*sizeof(float));
			if(allcandidates==NULL) return(1);
			
			size_t current_pos = 0;
			for(size_t f=0; f<nSubLists; f++){
				memcpy(&allcandidates[current_pos],clp[f]->candidates,clp[f]->nCandidates*4*sizeof(float));
				current_pos = current_pos + clp[f]->nCandidates*4;
			}
			return(0);
		}
		
		int ExportToFile(char *filename){
			size_t nSubLists = clp.size();
			size_t nCandidates_total = 0;
			for(size_t f=0; f<nSubLists; f++){
				nCandidates_total = nCandidates_total + clp[f]->nCandidates;
			}
			
			if(nCandidates_total>0){
				FILE *fp_out;
				if (( fp_out = fopen(filename, "wb") ) == NULL) return(1);
				
				size_t current_pos = 0;
				for(size_t f=0; f<nSubLists; f++){
					if(clp[f]->nCandidates>0){
						fwrite(clp[f]->candidates, clp[f]->nCandidates*sizeof(float), 4, fp_out);
					}
					current_pos = current_pos + clp[f]->nCandidates*4;
				}
				fclose(fp_out);
			}
			
			return(nCandidates_total);
		}
		
		JERK_Candidate getCandidate(size_t index){
			JERK_Candidate temp;
			if(index<nCandidates){
				temp.z   = allcandidates[index*4 + 0];
				temp.r   = allcandidates[index*4 + 1];
				temp.snr = allcandidates[index*4 + 2];
				temp.w   = allcandidates[index*4 + 3];
				return(temp);
			}
			else {
				temp.z   = 0;
				temp.r   = 0;
				temp.snr = 0;
				temp.w   = 0;
			}
			return(temp);
		}
		
		int clearSubLists(){
			for(size_t f=0; f<clp.size(); f++){
				clp[f]->~JERK_CandidateSubList();
			}
			clp.clear();
			return(0);
		}
		
		int clearCandidates(){
			if(allcandidates!=NULL) free(allcandidates);
			nCandidates = 0;
			return(0);
		}
		
		JERK_CandidateList(int reserve){
			allcandidates = NULL;
			nCandidates = 0;
			clp.reserve(reserve);
		}
		
		~JERK_CandidateList(){
			clearSubLists();
			clearCandidates();
		}
	};
	
} // namespace astroaccelerate
#endif
