#include "AstroAccelerate/headers_mains.h"
#include "device_dedispersion_kernel.cu"
#include "device_corner_turn_kernel.cu"
#include "device_binning_kernel.cu"

#include "device_SPS_inplace_kernel.cu" //Added by KA
#include "device_MSD_grid_kernel.cu" //Added by KA
#include "device_MSD_plane_kernel.cu" //Added by KA
#include "device_MSD_limited_kernel.cu" //Added by KA
#include "device_SNR_limited_kernel.cu" //Added by KA
#include "device_threshold_kernel.cu" //Added by KA
#include "device_single_FIR_kernel.cu" //Added by KA
#include "device_SPS_inplace.cu" //Added by KA
#include "device_MSD_grid.cu" //Added by KA
#include "device_MSD_plane.cu" //Added by KA
#include "device_MSD_limited.cu" //Added by KA
#include "device_SNR_limited.cu" //Added by KA
#include "device_threshold.cu" //Added by KA
#include "device_single_FIR.cu" //Added by KA

#include "device_bin.cu"
#include "device_dedisperse.cu"
#include "device_corner_turn.cu"
#include "device_set_stretch.cu"
#include "device_stats.cu"
#include "device_stretch.cu"
#include "device_power.cu"
#include "device_init.cu"
#include "device_inference.cu"
#include "device_load_data.cu"
#include "device_save_data.cu"
#include "device_zero_dm.cu"
#include "device_analysis.cu" //Added by KA
