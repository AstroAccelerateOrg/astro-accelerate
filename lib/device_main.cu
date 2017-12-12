#include "headers/headers_mains.h"
#include "device_dedispersion_kernel.cu"
#include "device_corner_turn_kernel.cu"
#include "device_binning_kernel.cu"


//MSD
#include "device_MSD_shared_kernel_functions.cu"
#include "device_MSD_normal_kernel.cu"
#include "device_MSD_outlier_rejection_kernel.cu"
#include "device_MSD.cu"
#include "device_MSD_legacy_kernel.cu"
#include "device_MSD_legacy.cu"



#include "device_SPS_inplace_kernel.cu" //Added by KA
#include "device_SNR_limited_kernel.cu" //Added by KA
#include "device_threshold_kernel.cu" //Added by KA
#include "device_single_FIR_kernel.cu" //Added by KA
#include "device_harmonic_summing_kernel.cu" //Added by KA


#include "device_SPS_inplace.cu" //Added by KA
#include "device_SPS_long.cu" //Added by KA
#include "device_SNR_limited.cu" //Added by KA
#include "device_threshold.cu" //Added by KA
#include "device_single_FIR.cu" //Added by KA
#include "device_harmonic_summing.cu"

#include "device_peak_find.cu"
#include "device_peak_find_kernel.cu"

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
#include "device_zero_dm_outliers.cu"
#include "device_rfi.cu"
#include "device_analysis.cu" //Added by KA
#include "device_periods.cu" //Added by KA

// fdas
#include "device_acceleration_fdas.cu"
