#ifndef __ASTROACCELERATE_SPS_CANDIDATELIST__
#define __ASTROACCELERATE_SPS_CANDIDATELIST__

struct SPS_Candidate {
	float DM;
	float time; //!!
	float SNR;
	float width;
};

class SPS_CandidateSubList;
class SPS_CandidateList;

class SPS_CandidateSubList { // time-DM data type
friend class SPS_CandidateList;
private:
	int nCandidates;
	int nBoxcars;
	
public:
	// data
	float *candidates;
	float *stdev;
	int   *boxcar_widths;
	
	// Coordinate transformation [indexes] -> [time, DM]
	float time_start;
	float sampling_time;
	float dm_step;
	float dm_low;
	float dm_high;
	int inBin;
	
	int getNumberOfCandidates(void) { return(nCandidates); }
	int getNumberOfBoxcars(void) { return(nBoxcars); }
	
	int exportToFile(void){
		char filename[200];
		
		sprintf(filename, "sublist-t_%.2f-dm_%.2f-%.2f.dat", time_start, dm_low, dm_high);
		
		FILE *fp_out;
		
		if(nCandidates>0){
			if (( fp_out = fopen(filename, "wb") ) == NULL) return(1);
			fwrite(candidates, nCandidates*sizeof(float), 4, fp_out);
			fclose(fp_out);
		}
		return(0);
	}
	
	SPS_CandidateSubList(int t_nCandidates, int t_nBoxcars, float *candidates_input, float *MSD_input, float *boxcar_widths_input){
		nCandidates = t_nCandidates;
		nBoxcars    = t_nBoxcars;
		
		if(nCandidates>0) {
			candidates = (float *) malloc(nCandidates*4*sizeof(float));
			if(candidates!=NULL) 
				memcpy(candidates, candidates_input, nCandidates*4*sizeof(float));
		}
		if(nBoxcars>0){
			stdev = (float *) malloc(nBoxcars*2*sizeof(float));
			boxcar_widths = (int *) malloc(nBoxcars*sizeof(int));
			if(stdev!=NULL) 
				memcpy(stdev, MSD_input, nBoxcars*2*sizeof(float));
			if(boxcar_widths!=NULL) 
				memcpy(boxcar_widths, boxcar_widths_input, nBoxcars*sizeof(int));
		}
	}
	
	~SPS_CandidateSubList(){
		if(candidates!=NULL) free(candidates);
		if(stdev!=NULL) free(stdev);
		if(boxcar_widths!=NULL) free(boxcar_widths);
	}
};

typedef SPS_CandidateSubList* CandidateSubListPointer;

class SPS_CandidateList {
private:

	
public:
	std::vector<CandidateSubListPointer> clp;
	float *allcandidates;
	size_t nCandidates;

	size_t getNumberOfCandidates(){
		return(nCandidates);
	}
	
	size_t getNumberOfSubLists(){
		return(clp.size());
	}
	
	SPS_CandidateSubList* getSubList(size_t index){
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
	
	int exportToFile(char *filename){
		FILE *fp_out;
		
		if(nCandidates>0){
			if (( fp_out = fopen(filename, "wb") ) == NULL) return(1);
			fwrite(allcandidates, nCandidates*sizeof(float), 4, fp_out);
			fclose(fp_out);
		}
		return(0);
	}
	
	SPS_Candidate getCandidate(size_t index){
		SPS_Candidate temp;
		temp.DM    = allcandidates[index*4 + 0];
		temp.time  = allcandidates[index*4 + 1];
		temp.SNR   = allcandidates[index*4 + 2];
		temp.width = allcandidates[index*4 + 3];
		return(temp);
	}
	
	int clearSubLists(){
		for(size_t f=0; f<clp.size(); f++){
			clp[f]->~SPS_CandidateSubList();
		}
		clp.clear();
		return(0);
	}
	
	int clearCandidates(){
		if(allcandidates!=NULL) free(allcandidates);
		nCandidates = 0;
		return(0);
	}
	
	SPS_CandidateList(int reserve){
		allcandidates = NULL;
		nCandidates = 0;
		clp.reserve(reserve);
	}
	
	~SPS_CandidateList(){
		clearSubLists();
		clearCandidates();
	}
};

#endif
