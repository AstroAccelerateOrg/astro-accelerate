#include <stdio.h>
#include <stdlib.h>

#include "headers/params.h"
#include "headers/host_MSD_stratagy.h"
#include "headers/device_MSD_plane_profile.h"
#include "headers/device_BC_plan.h"



//struct MSD_info{
//	unsigned long int max_timesamples;
//};


void Create_list_of_boxcar_widths2(std::vector<int> *boxcar_widths, std::vector<int> *BC_widths){
        int DIT_value, DIT_factor, width;
        DIT_value = 1;
        DIT_factor = 2;
        width = 0;
        for(int f=0; f<(int) BC_widths->size(); f++){
                for(int b=0; b<BC_widths->operator[](f); b++){
                        width = width + DIT_value;
                        boxcar_widths->push_back(width);
                }
                DIT_value = DIT_value*DIT_factor;
        }
}



void stratagy_MSD(int nDMs, float max_boxcar_width_in_sec, float tsamp, int nTimesamples, unsigned long int *maxtimesamples, size_t *MSD_profile_size_in_bytes, int *MSD_nDIT_widths){

        size_t vals;
        // Calculate the total number of values
        vals = ((size_t) nDMs)*((size_t) nTimesamples*NUM_STREAMS);

        size_t free_mem,total_mem;

        cudaMemGetInfo(&free_mem,&total_mem);
        printf("  Memory required by boxcar filters:%0.3f MB\n",(4.5*vals*sizeof(float) + 2*vals*sizeof(ushort))/(1024.0*1024) );
        printf("  Memory available:%0.3f MB \n", ((float) free_mem)/(1024.0*1024.0) );
	*maxtimesamples = (free_mem*0.95)/(NUM_STREAMS*(5.5*sizeof(float) + 2*sizeof(ushort)));
	printf("  Max samples: :%lld\n", (*maxtimesamples));

        int DMs_per_cycle = (*maxtimesamples)/nTimesamples;
        int itemp = (int) (DMs_per_cycle/THR_WARPS_PER_BLOCK);
        DMs_per_cycle = itemp*THR_WARPS_PER_BLOCK;
	printf("\n\n\nDMs_per_cycle: %i",DMs_per_cycle);
	
        int t_BC_widths[10]={PD_MAXTAPS,16,16,16,8,8,8,8,8,8};
        std::vector<int> BC_widths(t_BC_widths,t_BC_widths+sizeof(t_BC_widths)/sizeof(int));
        std::vector<PulseDetection_plan> PD_plan;
        std::vector<int> h_boxcar_widths;
        Create_list_of_boxcar_widths2(&h_boxcar_widths, &BC_widths);

	size_t MSD_DIT_profile_size_in_bytes, workarea_size_in_bytes;
        Get_MSD_plane_profile_memory_requirements(MSD_profile_size_in_bytes, &MSD_DIT_profile_size_in_bytes, &workarea_size_in_bytes, nTimesamples*NUM_STREAMS, nDMs, &h_boxcar_widths);
	
	int variable;
	MSD_nDIT_widths = &variable;
//	*MSD_nDIT_widths = 20;
        *MSD_nDIT_widths = (int)(max_boxcar_width_in_sec/tsamp);

        printf("\tSize MSD: %zu \tSize workarea: %zu, int: %i\n",  *MSD_profile_size_in_bytes, workarea_size_in_bytes/1024/1024.0, *MSD_nDIT_widths);


//        float *d_MSD_interpolated = NULL;
//        float *d_MSD_DIT = NULL;
//        float *temporary_workarea;
//        cudaMalloc((void **) &d_MSD_interpolated, MSD_profile_size_in_bytes);
//        cudaMalloc((void **) &temporary_workarea, workarea_size_in_bytes);


//	size_t MSD_profile_size_in_bytes, MSD_DIT_profile_size_in_bytes, workarea_size_in_bytes;
//	Get_MSD_plane_profile_memory_requirements(&MSD_profile_size_in_bytes, &MSD_DIT_profile_size_in_bytes, &workarea_size_in_bytes, nTimesamples, nDMs, &h_boxcar_widths);


}
