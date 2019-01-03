#ifndef ASTRO_ACCELERATE_AA_COMPUTE_HPP
#define ASTRO_ACCELERATE_AA_COMPUTE_HPP

#include <set>
#include <string>

/**
 * \namespace astroaccelerate
 * \brief A project-wide namespace astroaccelerate is used to separate the implementation of the library from library users.
 */
namespace astroaccelerate {
  
  /**
   * \namespace astroaccelerate::aa_compute 
   * \brief \brief The namespace that contains module names and the concept of a pipeline consisting of modules.
   */
  namespace aa_compute {

    /** \enum debug
     * \brief Contains debug flags.
     */
    enum class debug : int {
      debug = 0,
	analysis
	};
    
    /** \enum modules
     * \brief Contains the selectable modules.
     */
    enum class modules : int {
      empty = 0,
	dedispersion,
	analysis,
	periodicity,
	fdas,
	};

    /** \enum module_option
     * \brief Contains options for modules.
     */
    enum class module_option : int {
      empty = 0, //< The trivial module.
	zero_dm,
	zero_dm_with_outliers,
	rfi,
	old_rfi,
	sps_baseline_noise,
	output_dmt, //< Switches on output of ddtr to disk.
	output_ffdot_plan, //< Switches on output of ffdot_plan to disk.
	output_fdas_list, //< Switches on output of fdas_list to disk.
	candidate_algorithm, //< Enables/disables the candidate_algorithm
	fdas_custom_fft, //< Switches on output of custom_fft.
	fdas_inbin, //< Switches on inbin for fdas.
	fdas_norm //< Switches on norm for fdas.
	};

    /** \brief Function to convert module types into strings so that the user can query the pipeline. */
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
      case modules::fdas:
	return "fdas";
	break;
      case modules::periodicity:
	return "periodicity";
	break;
      default:
	return "ERROR: Module name not found";
	break;
      }
    }
    
    typedef std::set<aa_compute::modules> pipeline; /**< A pipeline is a collection of modules. */
    typedef std::set<aa_compute::module_option> pipeline_detail; /**< All pipeline options are contained in the pipeline_detail. */
    
  } // namespace aa_compute
} // namespace astroaccelerate


#endif // ASTRO_ACCELERATE_AA_COMPUTE_HPP
