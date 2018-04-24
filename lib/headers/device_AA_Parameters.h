#ifndef __ASTROACCELERATE_AA_PARAMETERS__
#define __ASTROACCELERATE_AA_PARAMETERS__

class AA_Parameters {
public:
	// Switches
	int enable_debug;
	int enable_analysis;
	int enable_acceleration;
	int enable_periodicity;
	int enable_zero_dm;
	int enable_zero_dm_with_outliers;
	int enable_rfi;
	int enable_old_rfi;
	int analysis_debug;
	int failsafe;
	
	AA_Parameters(int t_enable_debug, int t_enable_analysis, int t_enable_acceleration,	int t_enable_periodicity, int t_enable_zero_dm, int t_enable_zero_dm_with_outliers, int t_enable_rfi, int t_enable_old_rfi, int t_analysis_debug, int t_failsafe){
		enable_debug        = t_enable_debug;
		enable_analysis     = t_enable_analysis;
		enable_acceleration = t_enable_acceleration;
		enable_periodicity  = t_enable_periodicity;
		enable_zero_dm      = t_enable_zero_dm;
		enable_zero_dm_with_outliers = t_enable_zero_dm_with_outliers;
		enable_rfi          = t_enable_rfi;
		enable_old_rfi      = t_enable_old_rfi;
		analysis_debug      = t_analysis_debug;
		failsafe            = t_failsafe;
	}
};

#endif




	