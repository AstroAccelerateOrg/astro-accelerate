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
	int verbose;
	int export_data;
	int output_dedispersion_data;
	
	AA_Parameters(){
		enable_debug = 0;
		enable_analysis = 0;
		enable_acceleration = 0;
		enable_periodicity = 0;
		enable_zero_dm = 0;
		enable_zero_dm_with_outliers = 0;
		enable_rfi = 0;
		enable_old_rfi = 0;
		analysis_debug = 0;
		failsafe = 0;
		verbose = 1;
		export_data = 0;
		output_dedispersion_data = 0;
	}
	
	void Set(int t_enable_debug, int t_enable_analysis, int t_enable_acceleration,	int t_enable_periodicity, int t_enable_zero_dm, int t_enable_zero_dm_with_outliers, int t_enable_rfi, int t_enable_old_rfi, int t_analysis_debug, int t_failsafe){
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




	