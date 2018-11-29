#ifndef ASTRO_ACCELERATE_AA_COMPUTE_HPP
#define ASTRO_ACCELERATE_AA_COMPUTE_HPP

#include <set>
#include <string>
namespace astroaccelerate {
  namespace aa_compute {
    enum class debug : int {
			    debug = 0,
			    analysis
    };
        
    enum class modules : int {
			      empty = 0,
			      dedispersion,
			      analysis,
			      periodicity,
			      acceleration,
			      fdas_custom_fft,
			      fdas_inbin,
			      fdas_norm
    };

    enum class module_option : int {
				    empty = 0,
				    zero_dm,
				    zero_dm_with_outliers,
				    rfi,
				    old_rfi,
				    sps_baseline_noise,
				    output_dmt,
				    output_ffdot_plan,
				    output_fdas_list,
				    candidate_algorithm
    };

    //Function to convert module types into strings so that the user can query the pipeline
    inline const std::string module_name(const aa_compute::modules &module) {
      switch (module) {
      case modules::empty:
	return "empty";
	break;
      case modules::dedispersion:
	return "dedispersion";
	break;
      case modules::analysis:
	return "analysis";
	break;
      case modules::acceleration:
	return "acceleration";
	break;
      case modules::periodicity:
	return "periodicity";
	break;
      case modules::fdas_custom_fft:
	return "fdas_custom_fft";
	break;
      case modules::fdas_inbin:
	return "fdas_inbin";
	break;
      case modules::fdas_norm:
	return "fdas_norm";
	break;
      default:
	return "ERROR: Module name not found";
	break;
      }
    }
    
    typedef std::set<aa_compute::modules> pipeline;
    typedef std::set<aa_compute::module_option> pipeline_detail;
    
  }//namespace aa_compute
}//namespace astroaccelerate


#endif // ASTRO_ACCELERATE_AA_COMPUTE_HPP
