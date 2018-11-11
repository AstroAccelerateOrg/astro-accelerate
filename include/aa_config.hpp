#ifndef ASTRO_ACCELERATE_AA_CONFIG_HPP
#define ASTRO_ACCELERATE_AA_CONFIG_HPP

#include <stdio.h>
#include <string.h>
#include <set>
#include <vector>
#include <wordexp.h>
#include <algorithm>

#include "aa_ddtr_plan.hpp"
#include "aa_compute.hpp"
#include "aa_permitted_pipelines.hpp"

namespace astroaccelerate {

struct aa_CLI {
    std::vector<std::string> input;
};

struct aa_DM {
    float low;
    float high;
};

struct aa_config_flags {
    float narrow;
    float wide;
    float aggression;
    float nsearch;
    float power;
    float sigma_cutoff;
    float sigma_constant;
    float max_boxcar_width_in_sec;
    float periodicity_sigma_cutoff; //Should this be int or float?
    
    int multi_file;
    int output_dmt;
    int nboots;
    int ntrial_bins;
    int navdms;
    int range;
    int candidate_algorithm;
    int nb_selected_dm;
    int failsafe;
    int periodicity_nHarmonics;
    
    std::vector<int> inBin;
    std::vector<int> outBin;
    std::vector<aa_DM> selected_dm;
    struct User_DM {
        float low;
        float high;
        float step;
    };
    std::vector<User_DM> user_dm;
    std::vector<aa_compute::modules> user_modules;
    std::vector<aa_compute::debug> user_debug;
};

class aa_config {
public:
    
    aa_config(const std::string &config_file_path, const aa_CLI &cli_input) : configure_from_file(true), fpath(config_file_path), user_cli(cli_input) {
        
    }
    
    aa_config(aa_compute::pipeline &user_pipeline) : configure_from_file(false), pipeline(std::move(user_pipeline)) {
        
    }
    
    bool set_sigma(const float &sigma_constant, const float &sigma_cutoff) {
        flg.sigma_constant  = sigma_constant;
        flg.sigma_cutoff    = sigma_cutoff;
        return true;
    }
    
  bool add_dispersion_measure(const float &low, const float &high, const float &step, const int &inBin, const int &outBin) {
        aa_ddtr_plan::dm tmp = {low, high, step, inBin, outBin};
        m_dm_ranges.push_back(tmp);
        ++(flg.range);
        return true;
    }
    
    //This specialisation of the setup method happens only if reading from an input_file
    const aa_compute::pipeline setup() {
        if(configure_from_file) {
            if(get_user_input(user_cli)) {
                if(aa_permitted_pipelines::is_permitted(pipeline)) {
                    return pipeline;
                }
                else {
                    //User input was read successfully, but pipeline is not permitted
                    const aa_compute::pipeline empty = {aa_compute::modules::empty};
                    return empty;
                }
            }
            else {
                //Problem reading input file
                const aa_compute::pipeline empty = {aa_compute::modules::empty};
                return empty;
            }
        }
        else {
            //Configure from pre-supplied pipeline
            if(aa_permitted_pipelines::is_permitted(pipeline)) {
                return pipeline;
            }
            else {
                //The pre-supplied pipeline is not valid
                const aa_compute::pipeline empty = {aa_compute::modules::empty};
                return empty;
            }
        }
    }
    
protected:
    bool configure_from_file;
    std::string fpath;
    aa_config_flags flg;  // configuration flags
    aa_compute::pipeline pipeline;
    aa_CLI user_cli;
    std::vector<aa_ddtr_plan::dm> m_dm_ranges;
    
    
    bool get_user_input(const aa_CLI &cli_input) {
        const size_t argc = user_cli.input.size();
        
        FILE *fp_in = NULL;   // Path of the input file (configuration)
        FILE *fp = NULL;      // Path to the input data file (fil file)
        
        char string[100];
        
        if (argc < 2)
        {
            fprintf(stderr, "Need input file.\n");
            return false;
        }
        else if (argc == 2 && strcmp(user_cli.input[1].c_str(), "-help") != 0)
        {
            if (( fp_in = fopen(user_cli.input[1].c_str(), "r") ) == NULL)
            {
                fprintf(stderr, "Invalid input file!\n");
                return false;
            }
            flg.range = 0;
            while (!feof(fp_in))
            {
                if ( fscanf(fp_in, "%s", string) == 0 )
                {
                    fprintf(stderr, "failed to read input\n");
                    return false;
                }
                if (strcmp(string, "range") == 0) ++flg.range;
                if (strcmp(string, "selected_dm") == 0) ++flg.nb_selected_dm;
            }
            rewind(fp_in);
            
            flg.user_dm.resize(flg.range);
            flg.outBin.resize(flg.range);
            flg.inBin.resize(flg.range);
            
            
            // temporary variables to read dm range
            float temp_low  = 0;
            float temp_high = 0;
            float temp_step = 0;
            int temp_in_bin = 0;
            int temp_out_bin= 0;
            
            // read dm range if enabled
            int i=0;
            while (!feof(fp_in))
            {
                if ( fscanf(fp_in, "%s %f %f %f %d %d\n", string, &temp_low, &temp_high, &temp_step, &temp_in_bin, &temp_out_bin) == 0 )
                {
                    fprintf(stderr, "failed to read input\n");
                    exit(0);
                }
                if (strcmp(string, "range") == 0)
                {
                    flg.user_dm[i] = {temp_low, temp_high, temp_step};
                    flg.inBin[i]  = temp_in_bin;
                    flg.outBin[i] = temp_out_bin;
                    ++i;
                }
            }
            rewind(fp_in);
            
            flg.selected_dm.reserve(flg.nb_selected_dm);
            aa_DM tmp_selected;
            while (!feof(fp_in))
            {
                if ( fscanf(fp_in, "%s %f %f\n", string, &tmp_selected.low, &tmp_selected.high) == 0 )
                {
                    fprintf(stderr, "failed to read input\n");
                    return false;
                }
                if (strcmp(string, "selected_dm") == 0)
                {
                    flg.selected_dm.push_back(tmp_selected);
                }
            }
            rewind(fp_in);
            
            while (!feof(fp_in))
            {
                if ( fscanf(fp_in, "%s", string) == 0 )
                {
                    fprintf(stderr, "failed to read input\n");
                    return false;
                }
                if (strcmp(string, "debug") == 0)
                    flg.user_debug.push_back(aa_compute::debug::debug);
                if (strcmp(string, "analysis") == 0)
                    flg.user_modules.push_back(aa_compute::modules::analysis);
                if (strcmp(string, "periodicity") == 0)
                    flg.user_modules.push_back(aa_compute::modules::periodicity);
                if (strcmp(string, "acceleration") == 0)
                    flg.user_modules.push_back(aa_compute::modules::acceleration);
                if (strcmp(string, "output_ffdot_plan") == 0)
                    flg.user_modules.push_back(aa_compute::modules::output_ffdot_plan);
                if (strcmp(string, "output_fdas_list") == 0)
                    flg.user_modules.push_back(aa_compute::modules::output_fdas_list);
                if (strcmp(string, "output_dmt") == 0)
                    flg.user_modules.push_back(aa_compute::modules::output_dmt);
                if (strcmp(string, "zero_dm") == 0)
                    flg.user_modules.push_back(aa_compute::modules::zero_dm);
                if (strcmp(string, "zero_dm_with_outliers") == 0)
                    flg.user_modules.push_back(aa_compute::modules::zero_dm_with_outliers);
                if (strcmp(string, "rfi") == 0)
                    flg.user_modules.push_back(aa_compute::modules::rfi);
                if (strcmp(string, "oldrfi") == 0)
                    flg.user_modules.push_back(aa_compute::modules::old_rfi);
                if (strcmp(string, "threshold") == 0)
                    flg.user_modules.push_back(aa_compute::modules::candidate_algorithm);
                if (strcmp(string, "baselinenoise") == 0)
                    flg.user_modules.push_back(aa_compute::modules::sps_baseline_noise);
                if (strcmp(string, "fdas_custom_fft") == 0)
                    flg.user_modules.push_back(aa_compute::modules::fdas_custom_fft);
                if (strcmp(string, "fdas_inbin") == 0)
                    flg.user_modules.push_back(aa_compute::modules::fdas_inbin);
                if (strcmp(string, "fdas_norm") == 0)
                    flg.user_modules.push_back(aa_compute::modules::fdas_norm);
                if (strcmp(string, "multi_file") == 0)
                    flg.multi_file = 1;
                if (strcmp(string, "analysis_debug") == 0)
                    flg.user_debug.push_back(aa_compute::debug::analysis);
                if (strcmp(string, "failsafe") == 0)
                    flg.failsafe = 1;
                if (strcmp(string, "max_boxcar_width_in_sec") == 0)
                {
                    if ( fscanf(fp_in, "%f", &flg.max_boxcar_width_in_sec) == 0 )
                    {
                        fprintf(stderr, "failed to read input\n");
                        return false;
                    }
                }
                if (strcmp(string, "sigma_cutoff") == 0)
                {
                    if ( fscanf(fp_in, "%f", &flg.sigma_cutoff) == 0 )
                    {
                        fprintf(stderr, "failed to read input\n");
                        return false;
                    }
                }
                
                if (strcmp(string, "periodicity_sigma_cutoff") == 0) {
                    if ( fscanf(fp_in, "%f", &flg.periodicity_sigma_cutoff) == 0 ) {
                        fprintf(stderr, "failed to read input\n");
                        return false;
                    }
                }
                if (strcmp(string, "periodicity_harmonics") == 0) {
                    if ( fscanf(fp_in, "%d", &flg.periodicity_nHarmonics) == 0 ) {
                        fprintf(stderr, "failed to read input\n");
                        return false;
                    }
                }
                
                if (strcmp(string, "sigma_constant") == 0)
                {
                    if ( fscanf(fp_in, "%f", &flg.sigma_constant) == 0 )
                    {
                        fprintf(stderr, "failed to read input\n");
                        return false;
                    }
                }
                if (strcmp(string, "narrow") == 0)
                {
                    if ( fscanf(fp_in, "%f", &flg.narrow) == 0 )
                    {
                        fprintf(stderr, "failed to read input\n");
                        return false;
                    }
                }
                if (strcmp(string, "wide") == 0)
                {
                    if ( fscanf(fp_in, "%f", &flg.wide) == 0 )
                    {
                        fprintf(stderr, "failed to read input\n");
                        return false;
                    }
                }
                if (strcmp(string, "nboots") == 0)
                {
                    if ( fscanf(fp_in, "%d", &flg.nboots) == 0 )
                    {
                        fprintf(stderr, "failed to read input\n");
                        return false;
                    }
                }
                if (strcmp(string, "navdms") == 0)
                {
                    if ( fscanf(fp_in, "%d", &flg.navdms) == 0 )
                    {
                        fprintf(stderr, "failed to read input\n");
                        return false;
                    }
                }
                if (strcmp(string, "nwindows") == 0)
                {
                    if ( fscanf(fp_in, "%d", &flg.ntrial_bins) == 0 )
                    {
                        fprintf(stderr, "failed to read input\n");
                        return false;
                    }
                }
                if (strcmp(string, "aggression") == 0)
                {
                    if ( fscanf(fp_in, "%f", &flg.aggression) == 0 )
                    {
                        fprintf(stderr, "failed to read input\n");
                        return false;
                    }
                }
                if (strcmp(string, "power") == 0)
                {
                    if ( fscanf(fp_in, "%f", &flg.power) == 0 )
                    {
                        fprintf(stderr, "failed to read input\n");
                        return false;
                    }
                }
                if (strcmp(string, "file") == 0)
                {
                    // this command expand "~" to "home/username/"
                    wordexp_t expanded_string;
                    
                    if ( fscanf(fp_in, "%s", string) == 0 )
                    {
                        fprintf(stderr, "failed to read input\n");
                        return false;
                    }
                    wordexp(string, &expanded_string, 0);
                    if (( fp = fopen(expanded_string.we_wordv[0], "rb") ) == NULL)
                    {
                        fprintf(stderr, "Invalid data file!\n");
                        return false;
                    }
                    wordfree(&expanded_string);
                }
            }

	    //Apply overrides to the input txt file that were supplied via the CLI
	    if(std::find(cli_input.input.begin(), cli_input.input.end(), "narrow") != cli_input.input.end()) {
	      //Parse the override string...
	    }
	    
        }
        else if (argc == 2 && strcmp(user_cli.input[1].c_str(), "-help") == 0)
        {
            //help();
            return false;
        }
        else
        {
            fprintf(stderr, "Cannot recognise input, try \"./astro-accelerate -help.\"\n");
            return false;
        }
        
        return true;
    }
};

} //namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_CONFIG_HPP
