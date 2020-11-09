#ifndef ASTRO_ACCELERATE_AA_PIPELINE_API_HPP
#define ASTRO_ACCELERATE_AA_PIPELINE_API_HPP

#include <iostream>
#include <stdio.h>
#include <memory>
#include <map>

#include "aa_pipeline.hpp"
#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_analysis_plan.hpp"
#include "aa_analysis_strategy.hpp"
#include "aa_periodicity_plan.hpp"
#include "aa_periodicity_strategy.hpp"
#include "aa_jerk_plan.hpp"
#include "aa_jerk_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_device_info.hpp"
#include "aa_permitted_pipelines.hpp"
#include "aa_permitted_pipelines_1.hpp"
#include "aa_permitted_pipelines_2.hpp"
#include "aa_permitted_pipelines_3.hpp"
#include "aa_permitted_pipelines_3_0.hpp"
#include "aa_permitted_pipelines_4.hpp"
#include "aa_permitted_pipelines_4_0.hpp"
#include "aa_permitted_pipelines_5.hpp"
#include "aa_permitted_pipelines_5_0.hpp"
#include "aa_permitted_pipelines_generic.hpp"

#include "aa_log.hpp"

namespace astroaccelerate {

	/**
	 * \class aa_pipeline_api aa_pipeline_api.hpp "include/aa_pipeline_api.hpp"
	 * \brief Class to manage pipelines and their constituent components, plans, and strategies.
	 * \brief This class also delegates the movement of host memory into and out of the pipeline,
	 * \brief This class also manages the movement of device memory into and out of the device.
	 * \details The class is templated over the input and output data type.
	 * \details The class receives plan objects and calculates strategies.
	 * \details The user may obtain the strategies.
	 * \details The pipeline will not run unless all plans and strategies are successfully calculates for the pipeline that the user provided at construction.
	 * \details The pipeline strategy objects will request memory from the GPU that they will use when the pipeline is run.
	 * \details The construction of a new aa_pipeline_api object will reset all previously requested memory on the device (which may or may not have been allocated). \
	   Therefore, the next pipeline should be constructed after the previous one has been fully configured.
	 * \warning Configuring multiple pipeline objects and strategy objects at the same time means the pipeline will not see the correct amount of memory on the GPU.
	 * \todo Nice to have but not needed: Add a way to transfer ownership of the data between aa_pipeline_api objects.
	 * \author AstroAccelerate team
	 * \date: 23 October 2018.
	 */

	template<typename T>
	class aa_pipeline_api {
	private:
		std::map<aa_pipeline::component, bool> required_plans; /** Plans required to configure the pipeline. */
		std::map<aa_pipeline::component, bool> supplied_plans; /** Plans supplied by the user. */
		std::vector<aa_strategy*>              m_all_strategy; /** Base class pointers to all strategies bound to the pipeline. */
		aa_pipeline::pipeline                  m_requested_pipeline; /** The user requested pipeline that was bound to the aa_pipeline_api instance on construction. */
		const aa_pipeline::pipeline_option     m_pipeline_options; /** The user requested pipeline details containing component options for the aa_pipeline_api instance. */
		aa_device_info				           m_selected_device; /** The user provided GPU card information for the aa_pipeline_api instance. */
		std::unique_ptr<aa_pipeline_runner>    m_runner; /** A std::unique_ptr that will point to the correct class instantation of the selected aa_permitted_pipelines_ when the pipeline must be made ready to run. */

		aa_filterbank_metadata      m_filterbank_metadata; /** The filterbank file metadata that the user provided for the aa_pipeline_api instance on construction. */

		aa_ddtr_plan                m_ddtr_plan; /** The instance of this type that is currently bound to the aa_pipeline_api instance. */
		aa_ddtr_strategy            m_ddtr_strategy; /** The instance of this type that is currently bound to the aa_pipeline_api instance. */

		aa_analysis_plan            m_analysis_plan; /** The instance of this type that is currently bound to the aa_pipeline_api instance. */
		aa_analysis_strategy        m_analysis_strategy; /** The instance of this type that is currently bound to the aa_pipeline_api instance. */

		aa_periodicity_plan         m_periodicity_plan; /** The instance of this type that is currently bound to the aa_pipeline_api instance. */
		aa_periodicity_strategy     m_periodicity_strategy; /** The instance of this type that is currently bound to the aa_pipeline_api instance. */

		aa_fdas_plan                m_fdas_plan; /** The instance of this type that is currently bound to the aa_pipeline_api instance. */
		aa_fdas_strategy            m_fdas_strategy; /** The instance of this type that is currently bound to the aa_pipeline_api instance. */
		
		aa_jerk_plan                m_jerk_plan; /** The instance of this type that is currently bound to the aa_pipeline_api instance. */
		aa_jerk_strategy            m_jerk_strategy; /** The instance of this type that is currently bound to the aa_pipeline_api instance. */

		bool bound_with_raw_ptr; /** Flag to indicate whether the input data is bound via a raw pointer (true) or not (false). */
		bool pipeline_ready; /** Flag to indicate whether the pipeline is ready to execute (true) or not (false).  */

		std::vector<T>              data_in; /** Input data buffer. */
		T const*const               ptr_data_in; /** Input data pointer if bound_with_raw_ptr is true. */
		static int                  number_of_pipeline_instances;
		
	public:
	
		/** \brief Constructor for aa_pipeline_api that takes key parameters required on construction. */
		aa_pipeline_api(
			const aa_pipeline::pipeline &requested_pipeline,
			const aa_pipeline::pipeline_option &pipeline_options,
			const aa_filterbank_metadata &filterbank_metadata,
			T const*const input_data,
			aa_device_info &card_info)
			: 
			m_pipeline_options(pipeline_options),
			m_selected_device(card_info),
			m_filterbank_metadata(filterbank_metadata),
			bound_with_raw_ptr(true),
			pipeline_ready(false),
			ptr_data_in(input_data) 
		{
			//Add requested pipeline components
			for (auto i : requested_pipeline) {
				required_plans.insert(std::pair<aa_pipeline::component, bool>(i, true));
				supplied_plans.insert(std::pair<aa_pipeline::component, bool>(i, false));
			}
			m_all_strategy.reserve(requested_pipeline.size());
			m_requested_pipeline = requested_pipeline;

			++number_of_pipeline_instances;
			if (number_of_pipeline_instances > 1) {
				LOG(log_level::notice, "There is more than one aa_pipeline_api instance, make sure you constructed the next instance after the first instance is fully configured.");
				LOG(log_level::notice, "Otherwise, the new pipeline instance will not see the correct amount of GPU memory available.");
			}
		}

		/** \brief Destructor for aa_pipeline_api, performs cleanup as needed. */
		~aa_pipeline_api() {
			pipeline_ready = false;
			--number_of_pipeline_instances;
		}

		/**
		 * \brief Bind an aa_ddtr_plan to the aa_pipeline_api instance.
		 * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
		 */
		bool bind(aa_ddtr_plan plan) {
			pipeline_ready = false;

			//If a plan has already been supplied, return false and do nothing with the new plan
			if (supplied_plans.find(aa_pipeline::component::dedispersion) == supplied_plans.end()) {
				return false;
			}

			//Does the pipeline actually need this plan?
			if (required_plans.find(aa_pipeline::component::dedispersion) != required_plans.end()) {
				m_ddtr_plan = plan;

				//ddtr_strategy needs to know if analysis will be required
				if (required_plans.find(aa_pipeline::component::analysis) != required_plans.end()) {
					aa_ddtr_strategy ddtr_strategy(m_ddtr_plan, m_filterbank_metadata, m_selected_device.free_memory(), true, &m_selected_device);
					if (ddtr_strategy.ready()) {
						m_ddtr_strategy = std::move(ddtr_strategy);
						m_all_strategy.push_back(&m_ddtr_strategy);
						//If the plan is valid then the supplied_plan becomes true
						supplied_plans.at(aa_pipeline::component::dedispersion) = true;
					}
					else {
						return false;
					}
				}
				else {
					aa_ddtr_strategy ddtr_strategy(m_ddtr_plan, m_filterbank_metadata, m_selected_device.free_memory(), false, &m_selected_device);
					if (ddtr_strategy.ready()) {
						m_ddtr_strategy = std::move(ddtr_strategy);
						m_all_strategy.push_back(&m_ddtr_strategy);
						//If the plan is valid then the supplied_plan becomes true
						supplied_plans.at(aa_pipeline::component::dedispersion) = true;
					}
					else {
						return false;
					}
				}
			}
			else {
				//The plan is not required, ignore.
				return false;
			}

			return true;
		}


		/**
		 * \brief Bind an aa_analysis_plan to the aa_pipeline_api instance.
		 * \details If your pipeline includes analysis, then you must bind a valid aa_analysis_plan.
		 * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
		 */
		bool bind(aa_analysis_plan plan) {
			pipeline_ready = false;

			//If a plan has already been supplied, return false and do nothing with the new plan
			if (supplied_plans.find(aa_pipeline::component::analysis) == supplied_plans.end()) {
				return false;
			}

			//Does the pipeline actually need this plan?
			if (required_plans.find(aa_pipeline::component::analysis) != required_plans.end()) {
				//Is the ddtr_strategy provided by this analysis_plan ready?
				if (!plan.ddtr_strategy().ready()) {
					//This ddtr_strategy is not ready, so ignore this analysis_plan.
					printf("Not ready the strategy\n");
					return false;
				}

				m_analysis_plan = plan;

				aa_analysis_strategy analysis_strategy(m_analysis_plan, &m_selected_device);
				if (analysis_strategy.ready()) {
					m_analysis_strategy = std::move(analysis_strategy);
					m_all_strategy.push_back(&m_analysis_strategy);
				}
				else {
					return false;
				}

				//If the plan is valid then the supplied_plan becomes true
				supplied_plans.at(aa_pipeline::component::analysis) = true;
			}
			else {
				//The plan is not required, ignore.
				return false;
			}

			return true;
		}

		/**
		 * \brief Bind an aa_periodicity_plan to the aa_pipeline_api instance.
		 * \details If your pipeline includes periodicity, then you must bind a valid aa_periodicity_plan.
		 * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
		 */
		bool bind(aa_periodicity_plan plan) {
			pipeline_ready = false;

			//If a plan has already been supplied, return false and do nothing with the new plan
			if (supplied_plans.find(aa_pipeline::component::periodicity) == supplied_plans.end()) {
				return false;
			}

			//Does the pipeline actually need this plan?
			if (required_plans.find(aa_pipeline::component::periodicity) != required_plans.end()) {
				m_periodicity_plan = plan;
				aa_periodicity_strategy periodicity_strategy(m_periodicity_plan);
				if (periodicity_strategy.ready()) {
					m_periodicity_strategy = std::move(periodicity_strategy);
					m_all_strategy.push_back(&m_periodicity_strategy);
				}
				else {
					return false;
				}

				//If the plan is valid then the supplied_plan becomes true
				supplied_plans.at(aa_pipeline::component::periodicity) = true;
			}
			else {
				//The plan is not required, ignore.
				return false;
			}

			return true;
		}

		/**
		 * \brief Bind an aa_fdas_plan to the aa_pipeline_api instance.
		 * \details If your pipeline includes fdas, then you must bind a valid aa_fdas_plan.
		 * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
		 */
		bool bind(aa_fdas_plan plan) {
			pipeline_ready = false;

			//If a plan has already been supplied, return false and do nothing with the new plan
			if (supplied_plans.find(aa_pipeline::component::fdas) == supplied_plans.end()) {
				return false;
			}

			//Does the pipeline actually need this plan?
			if (required_plans.find(aa_pipeline::component::fdas) != required_plans.end()) {
				m_fdas_plan = plan;
				aa_fdas_strategy fdas_strategy(m_fdas_plan);
				if (fdas_strategy.ready()) {
					m_fdas_strategy = std::move(fdas_strategy);
					m_all_strategy.push_back(&m_fdas_strategy);
				}
				else {
					return false;
				}

				//If the plan is valid then the supplied_plan becomes true 
				supplied_plans.at(aa_pipeline::component::fdas) = true;
			}
			else {
				//The plan is not required, ignore.
				return false;
			}

			return true;
		}
		
		bool bind(aa_jerk_plan plan) {
			pipeline_ready = false;

			//If a plan has already been supplied, return false and do nothing with the new plan
			if (supplied_plans.find(aa_pipeline::component::jerk) == supplied_plans.end()) {
				return false;
			}

			//Does the pipeline actually need this plan?
			if (required_plans.find(aa_pipeline::component::jerk) != required_plans.end()) {
				m_jerk_plan = plan;
				aa_jerk_strategy jerk_strategy(m_jerk_plan);
				if (jerk_strategy.ready()) {
					m_jerk_strategy = std::move(jerk_strategy);
					m_all_strategy.push_back(&m_jerk_strategy);
				}
				else {
					return false;
				}

				//If the plan is valid then the supplied_plan becomes true 
				supplied_plans.at(aa_pipeline::component::jerk) = true;
			}
			else {
				//The plan is not required, ignore.
				return false;
			}

			return true;
		}

		/** \returns The aa_ddtr_strategy instance bound to the pipeline instance, or a trivial instance if a valid aa_ddtr_strategy does not yet exist. */
		aa_ddtr_strategy ddtr_strategy() {
			//Does the pipeline actually need this strategy?
			if (required_plans.find(aa_pipeline::component::dedispersion) != required_plans.end()) {
				//It does need this strategy.
				//Is it already computed?
				if (m_ddtr_strategy.ready()) {
					//Return since it was already computed.
					return m_ddtr_strategy;
				}
				else {
					//ddtr_strategy was not yet computed, do it now.
					//ddtr_strategy needs to know if analysis will be required
					if (required_plans.find(aa_pipeline::component::analysis) != required_plans.end()) { //analysis will be required
						aa_ddtr_strategy ddtr_strategy(m_ddtr_plan, m_filterbank_metadata, m_selected_device.free_memory(), true, &m_selected_device);
						if (ddtr_strategy.ready()) {
							m_ddtr_strategy = std::move(ddtr_strategy);
							m_all_strategy.push_back(&m_ddtr_strategy);
						}
						else { //Tried to calculate ddtr_strategy with analysis enabled, but failed.
							aa_ddtr_strategy empty_strategy;
							return empty_strategy;
						}
					}
					else { //analysis will not be required
						aa_ddtr_strategy ddtr_strategy(m_ddtr_plan, m_filterbank_metadata, m_selected_device.free_memory(), false, &m_selected_device);
						if (ddtr_strategy.ready()) {
							m_ddtr_strategy = std::move(ddtr_strategy);
							m_all_strategy.push_back(&m_ddtr_strategy);
						}
						else { //Tried to calculate ddtr_strategy with analysis disabled, but failed.
							aa_ddtr_strategy empty_strategy;
							return empty_strategy;
						}
					}
				}

				return m_ddtr_strategy;
			}
			else {
				//The pipeline does not need this strategy
				aa_ddtr_strategy empty_strategy;
				return empty_strategy;
			}
		}

		size_t get_nRanges(){
			return m_ddtr_strategy.get_nRanges();
		}

		const int* get_ndms_array(){
			return m_ddtr_strategy.ndms_data();
		}

		int dm_low(const int range){
			return m_ddtr_strategy.dm(range).low;
		}

		int total_computed_samples(){
			int tprocessed = 0;
			for(size_t j = 0; j < m_ddtr_strategy.t_processed().at(0).size(); j++) {
				tprocessed += m_ddtr_strategy.t_processed()[0][j];
			}
			return tprocessed;
		}

		int samples_in_current(int i, int j){
			return m_ddtr_strategy.t_processed()[i][j];
		}

		/** \returns The aa_analysis_strategy instance bound to the pipeline instance, or a trivial instance if a valid aa_analysis_strategy does not yet exist. */
		aa_analysis_strategy analysis_strategy() {
			//Does the pipeline actually need this strategy? 
			if (required_plans.find(aa_pipeline::component::analysis) != required_plans.end()) {
				//It does need this strategy.                                                                                                                                                                                        
				//Is it already computed?
				if (m_analysis_strategy.ready()) { //Return since it was already computed.
					return m_analysis_strategy;
				}
				else {
					//analysis_strategy was not yet computed, do it now.
					aa_analysis_strategy analysis_strategy(m_analysis_plan, &m_selected_device);
					if (analysis_strategy.ready()) {
						m_analysis_strategy = std::move(analysis_strategy);
						m_all_strategy.push_back(&m_analysis_strategy);
					}
					else { //Tried to calculate analysis_strategy, but failed.
						aa_analysis_strategy empty_strategy;
						return empty_strategy;
					}
				}
			}
			else {
				//The pipeline does not need this strategy
				aa_analysis_strategy empty_strategy;
				return empty_strategy;
			}
			aa_analysis_strategy empty_strategy;
			return empty_strategy;
		}

		/** \returns The aa_periodicity_strategy instance bound to the pipeline instance, or a trivial instance if a valid aa_periodicity_strategy does not yet exist. */
		aa_periodicity_strategy periodicity_strategy() {
			//Does the pipeline actually need this strategy?
			if (required_plans.find(aa_pipeline::component::periodicity) != required_plans.end()) {
				//It does need this strategy.
				//Is it already computed?
				if (m_periodicity_strategy.ready()) { //Return since it was already computed.
					return m_periodicity_strategy;
				}
				else {
					//periodicity_strategy was not yet computed, do it now.
					aa_periodicity_strategy periodicity_strategy(m_periodicity_plan);
					if (periodicity_strategy.ready()) {
						m_periodicity_strategy = std::move(periodicity_strategy);
						m_all_strategy.push_back(&m_periodicity_strategy);
					}
					else { //Tried to calculate periodicity strategy, but failed.
						aa_periodicity_strategy empty_strategy;
						return empty_strategy;
					}
				}
			}
			else {
				//The pipeline does not need this strategy
				aa_periodicity_strategy empty_strategy;
				return empty_strategy;

			}
		}

		/** \returns The aa_fdas_strategy instance bound to the pipeline instance, or a trivial instance if a valid aa_fdas_strategy does not yet exist. */
		aa_fdas_strategy fdas_strategy() {
			//Does the pipeline actually need this strategy?
			if (required_plans.find(aa_pipeline::component::fdas) != required_plans.end()) {
				//It does need this strategy.
				//Is it already computed?
				if (m_fdas_strategy.ready()) { //Return since it was already computed.
					return m_fdas_strategy;
				}
				else {
					//fdas_strategy was not yet computed, do it now.
					aa_fdas_strategy fdas_strategy(m_fdas_plan);
					if (fdas_strategy.ready()) {
						m_fdas_strategy = std::move(fdas_strategy);
						m_all_strategy.push_back(&m_fdas_strategy);
					}
					else { //Tried to calculate fdas strategy, but failed.
						aa_fdas_strategy empty_strategy;
						return empty_strategy;
					}
				}
			}
			else {
				//The pipeline does not need this strategy
				aa_fdas_strategy empty_strategy;
				return empty_strategy;
			}
			aa_fdas_strategy empty_strategy;
			return empty_strategy;
		}
		
		/** \returns The aa_jerk_strategy instance bound to the pipeline instance, or a trivial instance if a valid aa_jerk_strategy does not yet exist. */
		aa_jerk_strategy jerk_strategy() {
			//Does the pipeline actually need this strategy?
			if (required_plans.find(aa_pipeline::component::jerk) != required_plans.end()) {
				//It does need this strategy.
				//Is it already computed?
				if (m_jerk_strategy.ready()) { //Return since it was already computed.
					return m_jerk_strategy;
				}
				else {
					//fdas_strategy was not yet computed, do it now.
					aa_jerk_strategy jerk_strategy(m_jerk_plan);
					if (jerk_strategy.ready()) {
						m_jerk_strategy = std::move(jerk_strategy);
						m_all_strategy.push_back(&m_jerk_strategy);
					}
					else { //Tried to calculate fdas strategy, but failed.
						aa_jerk_strategy empty_strategy;
						return empty_strategy;
					}
				}
			}
			else {
				//The pipeline does not need this strategy
				aa_jerk_strategy empty_strategy;
				return empty_strategy;
			}
			aa_jerk_strategy empty_strategy;
			return empty_strategy;
		}

		/** \returns The aa_filterbank_metadata instance bound to the pipeline instance. */
		aa_filterbank_metadata metadata() {
			return m_filterbank_metadata;
		}

		/**
		 * \brief Check whether all components and strategies are ready for the pipeline to be able to execute.
		 * \returns A boolean to indicate whether the pipeline is ready (true) or not (false).
		 * \details If false, the user should check whether all plans have been bound and whether all strategies resulting from those plans are valid.
		 */
		bool ready() {
			if (pipeline_ready) {
				return true;
			}

			pipeline_ready = false;
			if (!aa_permitted_pipelines::is_permitted(m_requested_pipeline)) {
				LOG(log_level::error, "The requested pipeline is not permitted and cannot be made ready.");
				return false;
			}

			// Check if all plans are supplied.
			// If not, then one or more strategies will not be ready.
			// If so, return false.
			for (auto const& i : supplied_plans) {
				if (i.second == false) {
					LOG(log_level::error, aa_pipeline::component_name(i.first) + " plan is not ok.");
					return false;
				}
			}

			// Check if all strategies are ready.
			// If not, then return false.
			bool all_strategies_ready = true;
			for (auto strategy : m_all_strategy) {
				if (strategy->setup()) {
					//Memory allocations and setup happened successfully
				}
				else {
					//Setup failed.
					LOG(log_level::error, strategy->name()+" setup failed.");
					all_strategies_ready = false;
				}
			}

			if (!all_strategies_ready) {
				return false;
			}
			
			// By now pipeline is checked that it has reasonable components, that it has all required plans and that strategies are calculated correctly

			// Start configuring the pipeline that will be run.
			constexpr aa_pipeline::component_option zero_dm = aa_pipeline::component_option::zero_dm;
			constexpr aa_pipeline::component_option zero_dm_with_outliers = aa_pipeline::component_option::zero_dm_with_outliers;

			if (m_pipeline_options.find(zero_dm) == m_pipeline_options.end() &&
				m_pipeline_options.find(zero_dm_with_outliers) == m_pipeline_options.end()) {
				LOG(log_level::notice, "Neither zero_dm nor zero_dm_with_outliers were selected. Selection OFF.");
			}

			if (m_pipeline_options.find(zero_dm) != m_pipeline_options.end() &&
				m_pipeline_options.find(zero_dm_with_outliers) != m_pipeline_options.end()) {
				LOG(log_level::error, "Both zero_dm and zero_dm_with_outliers cannot simultaneously be selected.");
				return false;
			}

			//If fdas acceleration is used, these flags must be set
			bool fdas_enable_custom_fft = false;
			bool fdas_enable_inbin = false;
			bool fdas_enable_norm = false;
			bool fdas_enable_output_ffdot_plan = false;
			bool fdas_enable_output_list = false;

			if (m_pipeline_options.find(aa_pipeline::component_option::fdas_custom_fft) != m_pipeline_options.end()) {
				fdas_enable_custom_fft = true;
			}

			if (m_pipeline_options.find(aa_pipeline::component_option::fdas_inbin) != m_pipeline_options.end()) {
				fdas_enable_inbin = true;
			}

			if (m_pipeline_options.find(aa_pipeline::component_option::fdas_norm) != m_pipeline_options.end()) {
				fdas_enable_norm = true;
			}

			if (m_pipeline_options.find(aa_pipeline::component_option::output_ffdot_plan) != m_pipeline_options.end()) {
				fdas_enable_output_ffdot_plan = true;
			}

			if (m_pipeline_options.find(aa_pipeline::component_option::output_fdas_list) != m_pipeline_options.end()) {
				fdas_enable_output_list = true;
			}

			/**
			 * The following code requires some explanation, especially for developers unfamiliar with modern C++.
			 * The purpose of the following block of code is to be a runtime dispatch. What this does is look up
			 * which version of the pipeline function matches what the user requested.
			 *
			 * All pipelines are contained inside classes of the form aa_permitted_pipelines_x, where x is a number.
			 * These pipeline classes are templated over two parameters, a component_option relating to which zero_dm
			 * version will be used, and another to indicate which candidate algorithm will be used.
			 * The zero_dm option is contained inside aa_pipeline::component_option::zero_dm
			 * and aa_pipeline::component_option::zero_dm_with_outliers.
			 * To make the code less verbose, constexpr are used to make create "abbreviations" for this syntax.
			 *
			 * All aa_permitted_pipeline_x classes are derived from a base class aa_pipeline_runner.
			 * A std::unique_ptr is used because this provides automatic storage (when the std::unique_ptr goes out of
			 * scope, the destructor is called for the object to which it points).
			 * A base class pointer of type aa_pipeline_runner is used so that it can point to any of the aa_permitted_pipelines_x.
			 * Runtime polymorphism ensures the methods for the selected aa_permitted_pipeline_x class are called.
			 *
			 * The code checks the m_requested_pipeline against the possible pipelines, and selects the one matching the user
			 * settings. The syntax of std::unique_ptr implies that the templated class and the template parameters must be specified
			 * both when specifying the type of the std::unique_ptr, and when the type is "new"'d.
			 *
			 * Lastly, because the type is being constructed, the parameters for the constructor must also be provided.
			 *
			 */
			 
			//Launching generic pipeline which would call components from a permitted pipeline
			m_runner = std::unique_ptr<aa_permitted_pipelines_generic>(new aa_permitted_pipelines_generic(
					m_requested_pipeline,
					m_pipeline_options, 
					m_ddtr_strategy, 
					m_analysis_strategy, 
					m_periodicity_strategy, 
					m_fdas_strategy, 
					m_jerk_strategy, 
					fdas_enable_custom_fft, 
					fdas_enable_inbin, 
					fdas_enable_norm, 
					fdas_enable_output_ffdot_plan, 
					fdas_enable_output_list, 
					ptr_data_in));
			bool is_pipeline_set_to_runner = true;

			/*
			constexpr aa_pipeline::component_option off = aa_pipeline::component_option::empty;
			constexpr aa_pipeline::component_option old_rfi = aa_pipeline::component_option::old_rfi;
			constexpr bool use_old_rfi = true;
			constexpr bool no_rfi = false;
			//Check which pipeline the user has requested (given by m_requested_pipeline) against the possible permitted pipelines.
			//Then, assign a new object of that type to the base class pointer.
			if (m_requested_pipeline == aa_permitted_pipelines::pipeline1) {
				if (m_pipeline_options.find(zero_dm) != m_pipeline_options.end()) {
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_1<zero_dm, use_old_rfi>>(new aa_permitted_pipelines_1<zero_dm, use_old_rfi>(m_ddtr_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 1 with zero_dm, and old_rfi");
					}
					else {
						//details contain zero_dm and do not contain old_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_1<zero_dm, no_rfi>>(new aa_permitted_pipelines_1<zero_dm, no_rfi>(m_ddtr_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 1 with zero_dm, and no rfi");
					}
				}
				else if (m_pipeline_options.find(zero_dm_with_outliers) != m_pipeline_options.end()) {
					//details contain zero_dm_with_outliers
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm_with_outliers and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_1<zero_dm_with_outliers, use_old_rfi>>(new aa_permitted_pipelines_1<zero_dm_with_outliers, use_old_rfi>(m_ddtr_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 1 with zero_dm_with_outliers, and old_rfi");
					}
					else {
						//details contain zero_dm_with_outliers and do not contain older_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_1<zero_dm_with_outliers, no_rfi>>(new aa_permitted_pipelines_1<zero_dm_with_outliers, no_rfi>(m_ddtr_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 1 with zero_dm_with_outliers, and no rfi");
					}
				}
				else {
					LOG(log_level::notice, "Neither zero_dm nor zero_dm_with_outliers were specified in the options list. Selection OFF.");
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						m_runner = std::unique_ptr<aa_permitted_pipelines_1<off, use_old_rfi>>(new aa_permitted_pipelines_1<off, use_old_rfi>(m_ddtr_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 1 without zero_dm, and old_rfi");
					}
					else {
						m_runner = std::unique_ptr<aa_permitted_pipelines_1<off, no_rfi>>(new aa_permitted_pipelines_1<off, no_rfi>(m_ddtr_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 1 without zero_dm, and no rfi");
					}
				}
			}
			else if (m_requested_pipeline == aa_permitted_pipelines::pipeline2) {
				if (m_pipeline_options.find(zero_dm) != m_pipeline_options.end()) {
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_2<zero_dm, use_old_rfi>>(new aa_permitted_pipelines_2<zero_dm, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 2 with zero_dm, and old_rfi");
					}
					else {
						//details contain zero_dm and do not contain old_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_2<zero_dm, no_rfi>>(new aa_permitted_pipelines_2<zero_dm, no_rfi>(m_ddtr_strategy, m_analysis_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 2 with zero_dm, and no rfi");
					}
				}
				else if (m_pipeline_options.find(zero_dm_with_outliers) != m_pipeline_options.end()) {
					//details contain zero_dm_with_outliers
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm_with_outliers and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_2<zero_dm_with_outliers, use_old_rfi>>(new aa_permitted_pipelines_2<zero_dm_with_outliers, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 2 with zero_dm_with_outliers, and old_rfi");
					}
					else {
						//details contain zero_dm_with_outliers and do not contain older_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_2<zero_dm_with_outliers, no_rfi>>(new aa_permitted_pipelines_2<zero_dm_with_outliers, no_rfi>(m_ddtr_strategy, m_analysis_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 2 with zero_dm_with_outliers, and no rfi");
					}
				}
				else {
					LOG(log_level::notice, "Neither zero_dm nor zero_dm_with_outliers were specified in the options list. Selection OFF.");
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						m_runner = std::unique_ptr<aa_permitted_pipelines_2<off, use_old_rfi>>(new aa_permitted_pipelines_2<off, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 2 without zero_dm, and old_rfi");
					}
					else {
						m_runner = std::unique_ptr<aa_permitted_pipelines_2<off, no_rfi>>(new aa_permitted_pipelines_2<off, no_rfi>(m_ddtr_strategy, m_analysis_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 2 without zero_dm, and no rfi");
					}
				}
			}
			else if (m_requested_pipeline == aa_permitted_pipelines::pipeline3) {
				if (m_pipeline_options.find(zero_dm) != m_pipeline_options.end()) {
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_3<zero_dm, use_old_rfi>>(new aa_permitted_pipelines_3<zero_dm, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 3 with zero_dm, and old_rfi");
					}
					else {
						//details contain zero_dm and do not contain old_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_3<zero_dm, no_rfi>>(new aa_permitted_pipelines_3<zero_dm, no_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 3 with zero_dm, and no rfi");
					}
				}
				else if (m_pipeline_options.find(zero_dm_with_outliers) != m_pipeline_options.end()) {
					//details contain zero_dm_with_outliers
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm_with_outliers and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_3<zero_dm_with_outliers, use_old_rfi>>(new aa_permitted_pipelines_3<zero_dm_with_outliers, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 3 with zero_dm_with_outliers, and old_rfi");
					}
					else {
						//details contain zero_dm_with_outliers and do not contain older_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_3<zero_dm_with_outliers, no_rfi>>(new aa_permitted_pipelines_3<zero_dm_with_outliers, no_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 3 with zero_dm_with_outliers, and no rfi");
					}
				}
				else {
					LOG(log_level::notice, "Neither zero_dm nor zero_dm_with_outliers were specified in the options list. Selection OFF.");
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						m_runner = std::unique_ptr<aa_permitted_pipelines_3<off, use_old_rfi>>(new aa_permitted_pipelines_3<off, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 3 without zero_dm, and old_rfi");
					}
					else {
						m_runner = std::unique_ptr<aa_permitted_pipelines_3<off, no_rfi>>(new aa_permitted_pipelines_3<off, no_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 3 without zero_dm, and no rfi");
					}
				}
			}
			else if (m_requested_pipeline == aa_permitted_pipelines::pipeline4) {
				if (m_pipeline_options.find(zero_dm) != m_pipeline_options.end()) {
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_4<zero_dm, use_old_rfi>>(new aa_permitted_pipelines_4<zero_dm, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 4 with zero_dm, and old_rfi");
					}
					else {
						//details contain zero_dm and do not contain old_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_4<zero_dm, no_rfi>>(new aa_permitted_pipelines_4<zero_dm, no_rfi>(m_ddtr_strategy, m_analysis_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 4 with zero_dm, and no rfi");
					}
				}
				else if (m_pipeline_options.find(zero_dm_with_outliers) != m_pipeline_options.end()) {
					//details contain zero_dm_with_outliers
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm_with_outliers and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_4<zero_dm_with_outliers, use_old_rfi>>(new aa_permitted_pipelines_4<zero_dm_with_outliers, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 4 with zero_dm_with_outliers, and old_rfi");
					}
					else {
						//details contain zero_dm_with_outliers and do not contain older_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_4<zero_dm_with_outliers, no_rfi>>(new aa_permitted_pipelines_4<zero_dm_with_outliers, no_rfi>(m_ddtr_strategy, m_analysis_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 4 with zero_dm_with_outliers, and no rfi")
					}
				}
				else {
					LOG(log_level::notice, "Neither zero_dm nor zero_dm_with_outliers were specified in the options list. Selection OFF.");
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						m_runner = std::unique_ptr<aa_permitted_pipelines_4<off, use_old_rfi>>(new aa_permitted_pipelines_4<off, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 4 without zero_dm, and old_rfi");
					}
					else {
						m_runner = std::unique_ptr<aa_permitted_pipelines_4<off, no_rfi>>(new aa_permitted_pipelines_4<off, no_rfi>(m_ddtr_strategy, m_analysis_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 4 without zero_dm, and no rfi");
					}
				}
			}
			else if (m_requested_pipeline == aa_permitted_pipelines::pipeline5) {
				if (m_pipeline_options.find(zero_dm) != m_pipeline_options.end()) {
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_5<zero_dm, use_old_rfi>>(new aa_permitted_pipelines_5<zero_dm, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 5 with zero_dm, and old_rfi");
					}
					else {
						//details contain zero_dm and do not contain old_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_5<zero_dm, no_rfi>>(new aa_permitted_pipelines_5<zero_dm, no_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 5 with zero_dm, and no rfi");
					}
				}
				else if (m_pipeline_options.find(zero_dm_with_outliers) != m_pipeline_options.end()) {
					//details contain zero_dm_with_outliers
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm_with_outliers and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_5<zero_dm_with_outliers, use_old_rfi>>(new aa_permitted_pipelines_5<zero_dm_with_outliers, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 5 with zero_dm_with_outliers, and old_rfi");
					}
					else {
						//details contain zero_dm_with_outliers and do not contain older_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_5<zero_dm_with_outliers, no_rfi>>(new aa_permitted_pipelines_5<zero_dm_with_outliers, no_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 5 with zero_dm_with_outliers, and no rfi");
					}
				}
				else {
					LOG(log_level::notice, "Neither zero_dm nor zero_dm_with_outliers were specified in the options list. Selection OFF.");
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						m_runner = std::unique_ptr<aa_permitted_pipelines_5<off, use_old_rfi>>(new aa_permitted_pipelines_5<off, use_old_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 5 without zero_dm, and old_rfi");
					}
					else {
						m_runner = std::unique_ptr<aa_permitted_pipelines_5<off, no_rfi>>(new aa_permitted_pipelines_5<off, no_rfi>(m_ddtr_strategy, m_analysis_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 5 without zero_dm, and no rfi");
					}
				}
			}
			else if (m_requested_pipeline == aa_permitted_pipelines::pipeline3_0) {
				if (m_pipeline_options.find(zero_dm) != m_pipeline_options.end()) {
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_3_0<zero_dm, use_old_rfi>>(new aa_permitted_pipelines_3_0<zero_dm, use_old_rfi>(m_ddtr_strategy, m_periodicity_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 3_0 (no analysis) with zero_dm, and old_rfi");
					}
					else {
						//details contain zero_dm and do not contain old_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_3_0<zero_dm, no_rfi>>(new aa_permitted_pipelines_3_0<zero_dm, no_rfi>(m_ddtr_strategy, m_periodicity_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 3_0 (no analysis) with zero_dm, and no rfi");
					}
				}
				else if (m_pipeline_options.find(zero_dm_with_outliers) != m_pipeline_options.end()) {
					//details contain zero_dm_with_outliers
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm_with_outliers and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_3_0<zero_dm_with_outliers, use_old_rfi>>(new aa_permitted_pipelines_3_0<zero_dm_with_outliers, use_old_rfi>(m_ddtr_strategy, m_periodicity_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 3_0 (no analysis) with zero_dm_with_outliers, and old_rfi");
					}
					else {
						//details contain zero_dm_with_outliers and do not contain older_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_3_0<zero_dm_with_outliers, no_rfi>>(new aa_permitted_pipelines_3_0<zero_dm_with_outliers, no_rfi>(m_ddtr_strategy, m_periodicity_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 3_0 (no analysis) with zero_dm_with_outliers, and no rfi");
					}
				}
				else {
					LOG(log_level::notice, "Neither zero_dm nor zero_dm_with_outliers were specified in the options list. Selection OFF.");
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						m_runner = std::unique_ptr<aa_permitted_pipelines_3_0<off, use_old_rfi>>(new aa_permitted_pipelines_3_0<off, use_old_rfi>(m_ddtr_strategy, m_periodicity_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 3_0 (no analysis) without zero_dm, and old_rfi");
					}
					else {
						m_runner = std::unique_ptr<aa_permitted_pipelines_3_0<off, no_rfi>>(new aa_permitted_pipelines_3_0<off, no_rfi>(m_ddtr_strategy, m_periodicity_strategy, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 3_0 (no analysis) without zero_dm, and no rfi");
					}
				}
			}
			else if (m_requested_pipeline == aa_permitted_pipelines::pipeline4_0) {
				if (m_pipeline_options.find(zero_dm) != m_pipeline_options.end()) {
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_4_0<zero_dm, use_old_rfi>>(new aa_permitted_pipelines_4_0<zero_dm, use_old_rfi>(m_ddtr_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 4_0 (no analysis) with zero_dm, and old_rfi");
					}
					else {
						//details contain zero_dm and do not contain old_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_4_0<zero_dm, no_rfi>>(new aa_permitted_pipelines_4_0<zero_dm, no_rfi>(m_ddtr_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 4_0 (no analysis) with zero_dm, and no rfi");
					}
				}
				else if (m_pipeline_options.find(zero_dm_with_outliers) != m_pipeline_options.end()) {
					//details contain zero_dm_with_outliers
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm_with_outliers and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_4_0<zero_dm_with_outliers, use_old_rfi>>(new aa_permitted_pipelines_4_0<zero_dm_with_outliers, use_old_rfi>(m_ddtr_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 4_0 (no analysis) with zero_dm_with_outliers, and old_rfi");
					}
					else {
						//details contain zero_dm_with_outliers and do not contain older_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_4_0<zero_dm_with_outliers, no_rfi>>(new aa_permitted_pipelines_4_0<zero_dm_with_outliers, no_rfi>(m_ddtr_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 4_0 (no analysis) with zero_dm_with_outliers, and no rfi");
					}
				}
				else {
					LOG(log_level::notice, "Neither zero_dm nor zero_dm_with_outliers were specified in the options list. Selection OFF.");
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						m_runner = std::unique_ptr<aa_permitted_pipelines_4_0<off, use_old_rfi>>(new aa_permitted_pipelines_4_0<off, use_old_rfi>(m_ddtr_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 4_0 (no analysis) without zero_dm, and old_rfi");
					}
					else {
						m_runner = std::unique_ptr<aa_permitted_pipelines_4_0<off, no_rfi>>(new aa_permitted_pipelines_4_0<off, no_rfi>(m_ddtr_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 4_0 (no analysis) without zero_dm, and no rfi");
					}
				}
			}
			else if (m_requested_pipeline == aa_permitted_pipelines::pipeline5_0) {
				if (m_pipeline_options.find(zero_dm) != m_pipeline_options.end()) {
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_5_0<zero_dm, use_old_rfi>>(new aa_permitted_pipelines_5_0<zero_dm, use_old_rfi>(m_ddtr_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 5_0 (no analysis) with zero_dm, and old_rfi");
					}
					else {
						//details contain zero_dm and do not contain old_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_5_0<zero_dm, no_rfi>>(new aa_permitted_pipelines_5_0<zero_dm, no_rfi>(m_ddtr_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 5_0 (no analysis) with zero_dm, and no rfi");
					}
				}
				else if (m_pipeline_options.find(zero_dm_with_outliers) != m_pipeline_options.end()) {
					//details contain zero_dm_with_outliers
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						//details contain zero_dm_with_outliers and old_rfi
						m_runner = std::unique_ptr<aa_permitted_pipelines_5_0<zero_dm_with_outliers, use_old_rfi>>(new aa_permitted_pipelines_5_0<zero_dm_with_outliers, use_old_rfi>(m_ddtr_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 5_0 (no analysis) with zero_dm_with_outliers, and old_rfi");
					}
					else {
						//details contain zero_dm_with_outliers and do not contain older_rfi, so old_rfi is false
						m_runner = std::unique_ptr<aa_permitted_pipelines_5_0<zero_dm_with_outliers, no_rfi>>(new aa_permitted_pipelines_5_0<zero_dm_with_outliers, no_rfi>(m_ddtr_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 5_0 (no analysis) with zero_dm_with_outliers, and no rfi");
					}
				}
				else {
					LOG(log_level::notice, "Neither zero_dm nor zero_dm_with_outliers were specified in the options list. Selection OFF.");
					if (m_pipeline_options.find(old_rfi) != m_pipeline_options.end()) {
						m_runner = std::unique_ptr<aa_permitted_pipelines_5_0<off, use_old_rfi>>(new aa_permitted_pipelines_5_0<off, use_old_rfi>(m_ddtr_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 5_0 (no analysis) without zero_dm, and old_rfi");
					}
					else {
						m_runner = std::unique_ptr<aa_permitted_pipelines_5_0<off, no_rfi>>(new aa_permitted_pipelines_5_0<off, no_rfi>(m_ddtr_strategy, m_periodicity_strategy, m_fdas_strategy, fdas_enable_custom_fft, fdas_enable_inbin, fdas_enable_norm, fdas_enable_output_ffdot_plan, fdas_enable_output_list, ptr_data_in));
						is_pipeline_set_to_runner = true;
						LOG(log_level::notice, "Selected Pipeline 5_0 (no analysis) without zero_dm, and no rfi");
					}
				}
			}
			else {
				//Pipeline 0
			}
			*/
			
			//Do any last checks on the plans as a whole
			if (is_pipeline_set_to_runner) {
				pipeline_ready = true;
			}

			LOG(log_level::notice, "---PIPELINE DIAGNOSTIC INFORMATION---");
			m_selected_device.print_card_info();

			aa_filterbank_metadata::print_info(m_filterbank_metadata);

			if (required_plans.find(aa_pipeline::component::dedispersion) != required_plans.end()) {
				aa_ddtr_strategy::print_info(m_ddtr_strategy);
			}

			if (required_plans.find(aa_pipeline::component::analysis) != required_plans.end()) {
				aa_analysis_strategy::print_info(m_analysis_strategy);
			}

			if (required_plans.find(aa_pipeline::component::periodicity) != required_plans.end()) {
				aa_periodicity_strategy::print_info(m_periodicity_strategy);
			}

			if (required_plans.find(aa_pipeline::component::fdas) != required_plans.end()) {
				aa_fdas_strategy::print_info(m_fdas_strategy);
			}
			
			if (required_plans.find(aa_pipeline::component::jerk) != required_plans.end()) {
				aa_jerk_strategy::print_info(m_jerk_strategy);
			}

			pipeline_ready = true;
			return true;
		}


		/** \brief Runs the pipeline end-to-end. */
		bool run() {
			/**
			 * This method to be overloaded with all possible combinations of
			 * data that the user may wish to extract from any pipeline.
			 * Any methods that are not supported are compile-time errors because
			 * the base class must provide a method for it.
			 */
			if (pipeline_ready && m_runner->setup()) {
				aa_pipeline_runner::status status_code;
				while (m_runner->next(status_code)) {
					LOG(log_level::notice, "Pipeline running over next chunk.");
				}
				
				if(status_code==aa_pipeline_runner::status::error){
					LOG(log_level::notice, "Pipeline cannot proceed due to error.");
					return false;
				}
				else return true;
			}
			else {
				LOG(log_level::error, "Pipeline could not start/resume because either pipeline is not ready or runner is not setup.");
				return false;
			}
		}

		/**
		 * \brief Runs the pipeline one step at a time and provides a status code.
		 * \details The user must call this method like a call-back.
		 * \details The pipeline is finished when the return value is "false".
		 **/
		bool run(aa_pipeline_runner::status &status_code) {
			/**
			 * This method to be overloaded with all possible combinations of
			 * data that the user may wish to extract from any pipeline.
			 * Any methods that are not supported are compile-time errors because
			 * the base class must provide a method for it.
			 */
			if (pipeline_ready && m_runner->setup()) {
				LOG(log_level::notice, "Pipeline running over next chunk.");
				bool return_value = m_runner->next(status_code);
				if(status_code==aa_pipeline_runner::status::error){
					LOG(log_level::notice, "Pipeline cannot proceed due to error.");
					return false;
				}
				else return (return_value);
			}
			else {
				LOG(log_level::error, "Pipeline could not start/resume because either pipeline is not ready or runner is not setup.");
				status_code = aa_pipeline_runner::status::finished;
				return false;
			}
		}


		float* h_SPD_snr(){
			return m_runner->h_SPD_snr();
		}

		unsigned int* h_SPD_dm(){
			return m_runner->h_SPD_dm();
		}	

		unsigned int* h_SPD_width(){
			return m_runner->h_SPD_width();
		}	

		unsigned int* h_SPD_ts(){
			return m_runner->h_SPD_ts();
		}	

		size_t SPD_nCandidates(){
			return m_runner->get_SPD_nCandidates();
		}

		int get_current_range(){
			return m_runner->get_current_range();
		}

		int get_current_tchunk(){
			return m_runner->get_current_tchunk();
		}

		long int get_current_inc(){
			return m_runner->get_current_inc();
		}

		bool cleanup(){
			return m_runner->cleanup();
		}

		float ***output_buffer(){
			/**
			 * \brief Return the output of the DDTR. 
			 */
			if (m_pipeline_options.find(aa_pipeline::component_option::copy_ddtr_data_to_host) != m_pipeline_options.end()) {
				return m_runner->output_buffer();
			}
			else {
				LOG(log_level::error, "Could not get data from DDTR. The data are not copied from GPU memory to host memory. Enable option copy_DDTR_data_to_host.");
				return m_runner->output_buffer();
			}
		}		

		/**
		 * \brief Function pass input/output data from one aa_pipeline_api instance to another.
		 * \warning Not yet implemented.
		 */
		bool handoff(aa_pipeline_api &next_pipeline) {
			/**
			 * Handoff control over the data to the next pipeline.
			 *
			 * The user must always retrieve the data from the last pipeline in the chain.
			 */
			return true;
		}
	};

	template<typename T> int aa_pipeline_api<T>::number_of_pipeline_instances = 0;
} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PIPELINE_API_HPP

