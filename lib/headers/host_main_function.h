#ifndef ASTROACCELERATE_MAIN_FUNCTION_H_
#define ASTROACCELERATE_MAIN_FUNCTION_H_

void main_function (
	float **h_SPS_candidatelist,
	size_t *nSPScandidates,
	float ***output_buffer, //---> AA Output data on the HOST
	unsigned short  *input_buffer, //---> AA Input data on the HOST
	DDTR_Plan *DDTR_plan,
	AA_Parameters *AA_params,
	MSD_Parameters *MSD_params,
	SPS_Parameters *SPS_params,
	PRS_Parameters *PRS_params,
	FDAS_Parameters *FDAS_params,
	clock_t start_time // Time measurements
	);
#endif
