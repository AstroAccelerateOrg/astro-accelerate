#ifndef ASTRO_ACCELERATE_AA_PIPELINE_HPP
#define ASTRO_ACCELERATE_AA_PIPELINE_HPP

#include <set>
#include <string>

/**
 * \namespace astroaccelerate
 * \brief A project-wide namespace astroaccelerate is used to separate the implementation of the library from library users.
 */
namespace astroaccelerate {
  
  /**
   * \namespace astroaccelerate::aa_pipeline
   * \brief \brief The namespace that contains component names and the concept of a pipeline consisting of components.
   */
  namespace aa_pipeline {

    /** \enum debug
     * \brief Contains debug flags.
     */
    enum class debug : int {
      debug = 0,
	analysis
	};
    
    /** \enum component
     * \brief Contains the selectable components.
     */
    enum class component : int {
				empty = 0,
				dedispersion,
				analysis,
				periodicity,
				fdas,
    };
    
    /** \enum component_option
     * \brief Contains options for components.
     */
    enum class component_option : int {
				       empty = 0, //< The trivial component.
				       zero_dm,
				       zero_dm_with_outliers,
				       old_rfi,
				       msd_baseline_noise,
				       output_dmt, //< Switches on output of ddtr to disk.
				       output_ffdot_plan, //< Switches on output of ffdot_plan to disk.
				       output_fdas_list, //< Switches on output of fdas_list to disk.
				       candidate_algorithm, //< Enables/disables the candidate_algorithm
				       fdas_custom_fft, //< Switches on output of custom_fft.
				       fdas_inbin, //< Switches on inbin for fdas.
				       fdas_norm //< Switches on norm for fdas.
    };
    
    /** \brief Function to convert component types into strings so that the user can query the pipeline. */
    inline const std::string component_name(const aa_pipeline::component &component) {
      switch (component) {
      case aa_pipeline::component::empty:
	return "empty";
	break;
      case aa_pipeline::component::dedispersion:
	return "dedispersion";
	break;
      case aa_pipeline::component::analysis:
	return "analysis";
	break;
      case aa_pipeline::component::fdas:
	return "fdas";
	break;
      case aa_pipeline::component::periodicity:
	return "periodicity";
	break;
      default:
	return "ERROR: Component name not found";
	break;
      }
    }
    
    typedef std::set<aa_pipeline::component> pipeline; /**< A pipeline is a collection of components. */
    typedef std::set<aa_pipeline::component_option> pipeline_option; /**< All pipeline options contained in pipeline_option. */
    
  } // namespace aa_pipeline
} // namespace astroaccelerate


#endif // ASTRO_ACCELERATE_AA_PIPELINE_HPP
