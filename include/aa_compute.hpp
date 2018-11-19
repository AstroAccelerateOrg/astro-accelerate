#ifndef ASTRO_ACCELERATE_AA_COMPUTE_HPP
#define ASTRO_ACCELERATE_AA_COMPUTE_HPP

#include <set>
#include <string>

namespace aa_compute {
    enum class debug : int {
        debug = 0,
        analysis
    };
        
    enum class modules : int {
        empty = 0,
        dedispersion,
        analysis,
        acceleration,
        periodicity,
        dmt,
        zero_dm,
        zero_dm_with_outliers,
        rfi,
        old_rfi,
        sps_baseline_noise,
        fdas_custom_fft,
        fdas_inbin,
        fdas_norm,
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
            case modules::dmt:
                return "dmt";
                break;
            case modules::zero_dm:
                return "zero_dm";
                break;
            case modules::zero_dm_with_outliers:
                return "zero_dm_with_outliers";
                break;
            case modules::rfi:
                return "rfi";
                break;
            case modules::old_rfi:
                return "old_rfi";
                break;
            case modules::sps_baseline_noise:
                return "sps_baseline_noise";
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
            case modules::output_ffdot_plan:
                return "output_ffdot_plan";
                break;
            case modules::output_fdas_list:
                return "output_fdas_list";
                break;
            case modules::candidate_algorithm:
                return "candidate_algorithm";
                break;
            default:
                return "ERROR: Module name not found";
                break;
        }
    }
    
    typedef std::set<aa_compute::modules> pipeline;
        
}


#endif // ASTRO_ACCELERATE_AA_COMPUTE_HPP
