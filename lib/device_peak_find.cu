// James Sharpe's peak finding code

#include "AstroAccelerate/device_BC_plan.h"

#include <npp.h>
#include "helper_cuda.h"

#include "AstroAccelerate/params.h"
#include "device_peak_find_kernel.cu"


void PEAK_FIND(float *d_output_SNR, ushort *d_output_taps, float *d_peak_list, int nDMs, int nTimesamples, float threshold, int max_peak_size, int *gmem_peak_pos, int shift, std::vector<PulseDetection_plan> *PD_plan, int max_iteration){
	int output_offset, decimated_timesamples, local_offset;
	
	
	dim3 blockDim(32, 2, 1);
	dim3 gridSize(1 + ((nTimesamples-1)/blockDim.x), 1 + ((nDMs-1)/blockDim.y), 1);

	
	output_offset=0;
	for(int f=0; f<max_iteration; f++){
		decimated_timesamples = (nTimesamples>>f);
		//local_offset = offset>>f;
		local_offset = PD_plan->operator[](f).unprocessed_samples;
		if( (decimated_timesamples-local_offset-1)>0 ){		
			gridSize.x = 1 + ((decimated_timesamples-local_offset-1)/blockDim.x);
			gridSize.y = 1 + ((nDMs-1)/blockDim.y);
			gridSize.z = 1;		
			
			dilate_peak_find<<<gridSize, blockDim>>>(&d_output_SNR[output_offset], &d_output_taps[output_offset], d_peak_list, decimated_timesamples, nDMs, local_offset, threshold, max_peak_size, gmem_peak_pos, shift, (1<<f));
			output_offset = output_offset + nDMs*decimated_timesamples;
		}
	}
}