#ifndef ASTRO_ACCELERATE_AA_CONFIG_HPP
#define ASTRO_ACCELERATE_AA_CONFIG_HPP

#include <stdio.h>
#include <string.h>
#include <set>
#include <vector>
#include <wordexp.h>
#include <algorithm>

#include "aa_log.hpp"
#include "aa_params.hpp"
#include "aa_ddtr_plan.hpp"
#include "aa_pipeline.hpp"
#include "aa_permitted_pipelines.hpp"
#include "aa_host_help.hpp"


namespace astroaccelerate {

	/** \struct aa_command_line_arguments
	 * \brief A struct to contain a std::vector of command line argument strings.
	 */
	struct aa_command_line_arguments {
		std::vector<std::string> input;
	};

	/** \struct aa_config_flags
	 * \brief A struct to contain all configurations from an input file.
	 */
	class aa_config_flags {
	public:
		float power;                    /**< Power setting. */
		float sigma_cutoff;             /**< Sigma threshold for candidate selection. */
		float sigma_constant;           /**< Sigma cutoff for outlier rejection. */
		float max_boxcar_width_in_sec;  /**< Boxcar width in seconds setting. */
		float periodicity_sigma_cutoff; /**< Periodicity sigma cutoff setting. Should this be int or float? */
		float z_max;  
		float z_step; 
		float w_max;  
		float w_step; 

		int multi_file;                 /**< Multi file setting. This looks like it is deprecated. */
		int output_dmt;                 /**< Enables or disables ddtr output to disk. */
		int nRanges;                      /**< Range setting incremented to be the total number of user selected dm ranges. */
		int candidate_algorithm;        /**< Enables or disables use of candidate algorithm for analysis. */
		int nb_selected_dm;             /**< Incremented to be the total number of user selected dm ranges. Looks like a legacy duplicate of range. */
		int failsafe;                   /**< Flag to select the failsafe algorithm for dedispersion. */
		int periodicity_nHarmonics;     /**< Number of harmonics setting for periodicity. */
		int selected_card_id;           /**< Selected card id on this machine. */
		int dered;                      /** Enable deredning. */
		bool rfi;                       /**< Enable (true) or disable (false) host RFI reduction of the input data. */
		std::vector<aa_pipeline::debug> user_debug; /**< std::vector of debug flags. */
		
		void init(){
			power = 0.0;
			sigma_cutoff = 0.0;
			sigma_constant = 0.0;
			max_boxcar_width_in_sec = 0.0;
			periodicity_sigma_cutoff = 0.0;
			
			z_max = 0.0;
			z_step = 0.0;
			w_max = 0.0;
			w_step = 0.0;
			
			multi_file = 0;
			output_dmt = 0;
			nRanges = 0;
			candidate_algorithm = 0;
			nb_selected_dm = 0;
			failsafe = 0;
			periodicity_nHarmonics = 0;
			selected_card_id = 0;
			dered = 0;
			rfi = false;
			user_debug = std::vector<aa_pipeline::debug>();
		}
		
		aa_config_flags(){
			init();
		}
	};

	/**
	 * \class aa_config aa_config.hpp "include/aa_config.hpp"
	 * \brief Class to set up configuration flags.
	 * \brief The configuration can be defined by a library user.
	 * \brief Or, it can be configured from an input file.
	 * \details A configuration is required in order to construct a pipeline object.
	 * \details The pipeline object is described in aa_pipeline.
	 * \author Cees Carels.
	 * \date 5 November 2018.
	 */

	class aa_config {
	protected:
		bool configure_from_file; //Boolean flag to indicate on construction whether an input file will be used for configuration.
		std::string fpath; //Path to the input data file.
		aa_pipeline::pipeline        m_pipeline; //The pipeline object that is configured by an instance of this class using an input file. 
		aa_pipeline::pipeline_option m_pipeline_options; //The pipeline settings configured by an instance of this class using an input file.
		aa_command_line_arguments user_cli; //The user supplied command line argument settings.
		aa_config_flags flg;  //Configuration flags 
		aa_ddtr_plan m_ddtr_plan; //The ddtr_plan that is configured by an instance of this class using an input file.

		/** \brief Reads an input text file.
		 * \details Parses the string content of the file, sets the user flags, pipeline components, and pipeline details.
		 * \warning This function always adds aa_pipeline::component::dedispersion to the m_pipeline because the stnadalone always performs dedispersion.
		 */
		bool get_user_input() {
			const size_t argc = user_cli.input.size();

			FILE *fp_in = NULL;   // Path of the input file (configuration)
			FILE *fp = NULL;      // Path to the input data file (fil file)

			char string[100];

			if (argc < 2) {
				fprintf(stderr, "Need input file.\n");
				help();
				return false;
			}
			else if (argc == 2 && strcmp(user_cli.input[1].c_str(), "-help") != 0) {
				if ((fp_in = fopen(user_cli.input[1].c_str(), "r")) == NULL) {
					fprintf(stderr, "Invalid input file!\n");
					return false;
				}
				flg.nRanges = 0;
				while (!feof(fp_in))
				{
					if (fscanf(fp_in, "%s", string) == 0)
					{
						fprintf(stderr, "failed to read input\n");
						return false;
					}
					if (strcmp(string, "range") == 0) ++flg.nRanges;
					if (strcmp(string, "selected_dm") == 0) ++flg.nb_selected_dm;
				}
				rewind(fp_in);

				// temporary variables to read dm range
				float temp_low = 0;
				float temp_high = 0;
				float temp_step = 0;
				int temp_in_bin = 0;
				int temp_out_bin = 0;

				// read dm range if enabled
				while (!feof(fp_in)) {
					if (fscanf(fp_in, "%s %f %f %f %d %d\n", string, &temp_low, &temp_high, &temp_step, &temp_in_bin, &temp_out_bin) == 0) {
						fprintf(stderr, "failed to read input\n");
						return false;
					}
					if (strcmp(string, "range") == 0) {
						m_ddtr_plan.add_dm(temp_low, temp_high, temp_step, temp_in_bin, temp_out_bin);
					}
				}
				rewind(fp_in);
				m_pipeline.insert(aa_pipeline::component::dedispersion); //Always add dedispersion to the pipeline
				while (!feof(fp_in)) {
					if (fscanf(fp_in, "%s", string) == 0) {
						fprintf(stderr, "failed to read input\n");
						return false;
					}
					if (strcmp(string, "debug") == 0)
						flg.user_debug.push_back(aa_pipeline::debug::debug);
					if (strcmp(string, "analysis") == 0)
						m_pipeline.insert(aa_pipeline::component::analysis);
					if (strcmp(string, "periodicity") == 0)
						m_pipeline.insert(aa_pipeline::component::periodicity);
					if (strcmp(string, "acceleration") == 0)
						m_pipeline.insert(aa_pipeline::component::fdas);
					if (strcmp(string, "acceleration_jerk") == 0){
						m_pipeline.insert(aa_pipeline::component::jerk);
					}
					if (strcmp(string, "output_ffdot_plan") == 0)
						m_pipeline_options.insert(aa_pipeline::component_option::output_ffdot_plan);
					if (strcmp(string, "output_fdas_list") == 0)
						m_pipeline_options.insert(aa_pipeline::component_option::output_fdas_list);
					if (strcmp(string, "output_dmt") == 0)
						m_pipeline_options.insert(aa_pipeline::component_option::output_dmt);
					if (strcmp(string, "dered") == 0)
						m_pipeline_options.insert(aa_pipeline::component_option::dered);
					if (strcmp(string, "zero_dm") == 0)
						m_pipeline_options.insert(aa_pipeline::component_option::zero_dm);
					if (strcmp(string, "zero_dm_with_outliers") == 0)
						m_pipeline_options.insert(aa_pipeline::component_option::zero_dm_with_outliers);
					if (strcmp(string, "input_DDTR_normalization") == 0)
						m_pipeline_options.insert(aa_pipeline::component_option::input_DDTR_normalization);
					if (strcmp(string, "output_DDTR_normalization") == 0)
						m_pipeline_options.insert(aa_pipeline::component_option::output_DDTR_normalization);
					if (strcmp(string, "set_bandpass_average") == 0)
						m_pipeline_options.insert(aa_pipeline::component_option::set_bandpass_average);
					if (strcmp(string, "rfi") == 0)
						flg.rfi = true;
					if (strcmp(string, "oldrfi") == 0)
						m_pipeline_options.insert(aa_pipeline::component_option::old_rfi);
					if (strcmp(string, "threshold") == 0) {
						m_pipeline_options.insert(aa_pipeline::component_option::candidate_algorithm);
						flg.candidate_algorithm = 1;
					}
					if (strcmp(string, "baselinenoise") == 0)
						m_pipeline_options.insert(aa_pipeline::component_option::msd_baseline_noise);
					if (strcmp(string, "fdas_custom_fft") == 0)
						m_pipeline_options.insert(aa_pipeline::component_option::fdas_custom_fft);
					if (strcmp(string, "fdas_inbin") == 0)
						m_pipeline_options.insert(aa_pipeline::component_option::fdas_inbin);
					if (strcmp(string, "fdas_norm") == 0)
						m_pipeline_options.insert(aa_pipeline::component_option::fdas_norm);
					if (strcmp(string, "multi_file") == 0)
						flg.multi_file = 1;
					if (strcmp(string, "analysis_debug") == 0)
						flg.user_debug.push_back(aa_pipeline::debug::analysis);
					if (strcmp(string, "failsafe") == 0)
						flg.failsafe = 1;
					if (strcmp(string, "max_boxcar_width_in_sec") == 0) {
						if (fscanf(fp_in, "%f", &flg.max_boxcar_width_in_sec) == 0)	{
							fprintf(stderr, "failed to read input\n");
							return false;
						}
					}
					if (strcmp(string, "sigma_cutoff") == 0) {
						if (fscanf(fp_in, "%f", &flg.sigma_cutoff) == 0) {
							fprintf(stderr, "failed to read input\n");
							return false;
						}
					}
					if (strcmp(string, "periodicity_sigma_cutoff") == 0) {
						if (fscanf(fp_in, "%f", &flg.periodicity_sigma_cutoff) == 0) {
							fprintf(stderr, "failed to read input\n");
							return false;
						}
					}
					if (strcmp(string, "periodicity_harmonics") == 0) {
						if (fscanf(fp_in, "%d", &flg.periodicity_nHarmonics) == 0) {
							fprintf(stderr, "failed to read input\n");
							return false;
						}
					}
					if (strcmp(string, "z_max") == 0) {
						if (fscanf(fp_in, "%f", &flg.z_max) == 0) {
							fprintf(stderr, "failed to read input\n");
							return false;
						}
					}
					if (strcmp(string, "z_step") == 0) {
						if (fscanf(fp_in, "%f", &flg.z_step) == 0) {
							fprintf(stderr, "failed to read input\n");
							return false;
						}
					}
					if (strcmp(string, "w_max") == 0) {
						if (fscanf(fp_in, "%f", &flg.w_max) == 0) {
							fprintf(stderr, "failed to read input\n");
							return false;
						}
					}
					if (strcmp(string, "w_step") == 0) {
						if (fscanf(fp_in, "%f", &flg.w_step) == 0) {
							fprintf(stderr, "failed to read input\n");
							return false;
						}
					}
					if (strcmp(string, "sigma_constant") == 0)
					{
						if (fscanf(fp_in, "%f", &flg.sigma_constant) == 0)
						{
							fprintf(stderr, "failed to read input\n");
							return false;
						}
					}
					if (strcmp(string, "power") == 0)
					{
						if (fscanf(fp_in, "%f", &flg.power) == 0)
						{
							fprintf(stderr, "failed to read input\n");
							return false;
						}
						m_ddtr_plan.set_power(flg.power);
					}
					if (strcmp(string, "selected_card_id") == 0)
					{
						if (fscanf(fp_in, "%d", &flg.selected_card_id) == 0)
						{
							fprintf(stderr, "failed to read input\n");
							return false;
						}
					}
					if (strcmp(string, "file") == 0) {
						// This command expands "~" to "home/username/"
						wordexp_t expanded_string;

						if (fscanf(fp_in, "%s", string) == 0)
						{
							fprintf(stderr, "failed to read input\n");
							return false;
						}
						wordexp(string, &expanded_string, 0);
						if ((fp = fopen(expanded_string.we_wordv[0], "rb")) == NULL)
						{
							fprintf(stderr, "Invalid data file!\n");
							help();
							return false;
						}
						else {
							//File is valid
							//Close it again for later re-opening
							fpath = expanded_string.we_wordv[0];
							fclose(fp);
						}
						wordfree(&expanded_string);
					}
				}
			}
			else if ((argc == 2 && (strcmp(user_cli.input[1].c_str(), "-help") == 0))
				||
				(argc == 2 && strcmp(user_cli.input[1].c_str(), "--help") == 0)) {
				help();
				return false;
			}
			else {
				fprintf(stderr, "Cannot recognise input, try \"./astro-accelerate -help.\"\n");
				return false;
			}

			return true;
		}
	
	public:
		/** \brief Constructor for aa_config to configure from an input file for the standalone executable.
		 * \brief The aa_command_line_arguments parameter cli_input is used to supply all command line interface arguments.
		 * \details This must at least contain the full path to an input file.
		 * \details If using AstroAccelerate as a library user, then use the other constructor of the aa_config class.
		 */
		aa_config(const aa_command_line_arguments &cli_input) {
			configure_from_file = true;
			user_cli = cli_input;
			flg.init();
			flg.power = 2.0; // The initialiser list is rather long, and if new members are added, the change in declaration order may introduce a bug. So, it is done explicitly in the body.
			flg.periodicity_nHarmonics = 32; // wrong
			flg.selected_card_id = CARD; // wrong
		}

		/** \brief Constructor for aa_config to configure via a user requested aa_pipeline::pipeline object.
		 * \details The aa_pipeline::pipeline object contains the user requested pipeline component to be run.
		 * \details Use this constructor when using AstroAccelerate as a library user.
		 */
		aa_config(aa_pipeline::pipeline &user_pipeline) {
			configure_from_file = false;
			m_pipeline = user_pipeline;
			flg.init();
			flg.power = 2.0; // The initialiser list is rather long, and if new members are added, the change in declaration order may introduce a bug. So, it is done explicitly in the body.
			flg.periodicity_nHarmonics = 32; // wrong
			flg.selected_card_id = CARD; // wrong
		}

		/** \brief Overloaded function that simplifies the setup method.
		 * \details A library user can check the validity of a pipeline by creating a new aa_config object as a library user, and calling this method.
		 * \details Alternatively, the library user can choose not to use aa_config at all.
		 * \details This method cannot be used in conjunction with reading an input file for the stnadalone.
		 * \returns A pipeline object that is either valid or empty if it is not valid.
		 */
		const aa_pipeline::pipeline setup() {
			return setup(m_ddtr_plan, flg, m_pipeline_options, fpath);
		}

		/** \brief Sets up the pipeline flags and ddtr_plan.
		 * \details This method is useful mainly for the standalone code, since library users can configure the supplied parameters themselves.
		 * \details If a library user wishes to use an input file for configuration, then they use this setup method.
		 * \details The method checks whether the aa_config was constructed from an input file or not, and if so it sets the provided plan and flags.
		 * \details If the aa_config object was constructed without an input file path, then this module only checks whether the supplied pipeline is valid.
		 * \returns A pipeline object that is either valid or empty if it is not valid.
		 */
		const aa_pipeline::pipeline setup(aa_ddtr_plan &ddtr_plan, aa_config_flags &user_flags, aa_pipeline::pipeline_option &pipeline_options, std::string &file_path) {
			if (configure_from_file) {
				if (get_user_input()) { //Read user input from text input file and set flag object.
					if (aa_permitted_pipelines::is_permitted(m_pipeline)) { //get_user_input has configured m_pipeline, so now its validity can be checked.
						ddtr_plan = m_ddtr_plan;
						file_path = fpath;
						user_flags = flg;
						pipeline_options = m_pipeline_options;
						return m_pipeline;
					}
					else {
						//User input was read successfully, but pipeline is not permitted
						LOG(log_level::error, "Pipeline configuration is unknown.");
						const aa_pipeline::pipeline empty = { aa_pipeline::component::empty };
						return empty;
					}
				}
				else {
					//Problem reading input file
					const aa_pipeline::pipeline empty = { aa_pipeline::component::empty };
					return empty;
				}
			}
			else {
				//Configure from pre-supplied pipeline
				if (aa_permitted_pipelines::is_permitted(m_pipeline)) {
					return m_pipeline;
				}
				else {
					//The pre-supplied pipeline is not valid
					LOG(log_level::error, "Pipeline configuration is unknown.");
					const aa_pipeline::pipeline empty = { aa_pipeline::component::empty };
					return empty;
				}
			}
		}
		
	};

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_CONFIG_HPP
