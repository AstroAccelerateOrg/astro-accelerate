#ifndef ASTRO_ACCELERATE_AA_DEVICE_INFO_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_INFO_HPP

#include <cuda.h>
#include <cuda_runtime.h>

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#include "aa_log.hpp"

/**
 * If the build system (CMake) defines aa_version.hpp, then include it.
 */
#ifdef ASTRO_ACCELERATE_VERSION_H_DEFINED
#include "aa_version.hpp"
#endif

namespace astroaccelerate {

/**
 * \class aa_device_info aa_device_info.hpp "include/aa_device_info.hpp"
 * \brief Obtain information about available GPUs and select the GPU to use for data processing.
 * \author AstroAccelerate
 * \date 1 November 2018.
 */

class aa_device_info {
public:

	//aa_device_info(const aa_device_info &) = delete; /** \brief Delete copy and move constructors. */
	//aa_device_info(aa_device_info &&) = delete; /** \brief Delete copy and move constructors. */
	//aa_device_info &operator=(const aa_device_info &) = delete; /** \brief Delete copy and move constructors. */
	//aa_device_info &operator=(aa_device_info &&) = delete; /** \brief Delete copy and move constructors. */

	/**
	 * \struct aa_card_info
	 * \brief Struct to contain CUDA card information.
	 */
	struct aa_card_info {
		int compute_capability_major;
		int compute_capability_minor;
		float clock_rate; // [kHz]
		float memory_clock_rate; // [kHz]
		float memory_bus_width; // [bits]
		int card_number;
		int driver_version;
		int runtime_version;
		int l2_cache_size; // [bytes]
		size_t global_memory; // [bytes]
		size_t total_memory; // [bytes]
		size_t free_memory; // [bytes]
		size_t available_memory; // [bytes]
		size_t user_requested_memory_for_allocation; // Memory the user wishes to allocate, may or may not already have been allocated.
		std::string name;
	};

	aa_device_info(int selected_device, size_t maximum_allowed_memory_for_pipeline = 0) {
		is_init = false;
		user_requested_memory_for_allocation = 0;
		selected_device_id = selected_device;
		check_for_devices();
		init_card();
		if (is_init==true) {
			if (maximum_allowed_memory_for_pipeline == 0)  pipeline_memory_limit = selected_card_info.global_memory;
			else pipeline_memory_limit = maximum_allowed_memory_for_pipeline;
		}
	}

	/** \brief The destructor is public so static references returned by the instance() method are deleted. */
	~aa_device_info() {
		cudaDeviceReset();
	}

	/** \returns The currently available free memory on the currently selected GPU as reported by the CUDA driver. */
	size_t free_memory() {
		if (is_init) {
			size_t free = 0;
			size_t total = 0;
			cudaMemGetInfo(&free, &total);
			selected_card_info.free_memory = free;
			if (free >= pipeline_memory_limit) selected_card_info.available_memory = pipeline_memory_limit - user_requested_memory_for_allocation;
			else selected_card_info.available_memory = free - user_requested_memory_for_allocation;
			return (selected_card_info.available_memory);
		}

		return 0;
	}
	
	bool request_memory(size_t memory_size_in_bytes){
		if( selected_card_info.available_memory > memory_size_in_bytes ) {
			user_requested_memory_for_allocation += memory_size_in_bytes;
			selected_card_info.available_memory = selected_card_info.available_memory - memory_size_in_bytes;
			return true;
		}
		else return false;
	}
	
	bool requested_memory(){
		return user_requested_memory_for_allocation;
	}
	
	bool reset_requested_memory(){
		user_requested_memory_for_allocation = 0;
		return true;
	}

	int currently_selected_card_id() const {
		return selected_device_id;
	}

	/** \brief Static method for printing member data for an instance of aa_card_info. */
	bool print_card_info() {
		LOG(log_level::dev_debug, "CARD INFORMATION:");
		LOG(log_level::dev_debug, "Name:\t\t\t"+selected_card_info.name);
		LOG(log_level::dev_debug, "Compute capability:\t" + std::to_string(selected_card_info.compute_capability_major) + "." + std::to_string(selected_card_info.compute_capability_minor));
		LOG(log_level::dev_debug, "Clock rate [kHz]:\t" + std::to_string(selected_card_info.clock_rate));
		LOG(log_level::dev_debug, "Memory clock rate [kHz]:\t" + std::to_string(selected_card_info.memory_clock_rate));
		LOG(log_level::dev_debug, "Memory bus width [bits]:\t" + std::to_string(selected_card_info.memory_bus_width));
		LOG(log_level::dev_debug, "ID # on this machine:\t" + std::to_string(selected_card_info.card_number));
		LOG(log_level::dev_debug, "Driver version:\t\t" + std::to_string(selected_card_info.driver_version));
		LOG(log_level::dev_debug, "l2_cache_size [bytes]:\t" + std::to_string(selected_card_info.l2_cache_size));
		LOG(log_level::dev_debug, "Global memory [bytes]:\t" + std::to_string(selected_card_info.global_memory));
		LOG(log_level::dev_debug, "Total memory [bytes]:\t" + std::to_string(selected_card_info.total_memory));
		LOG(log_level::dev_debug, "Free memory [bytes]:\t" + std::to_string(selected_card_info.free_memory));
		return true;
	}

private:
	int selected_device_id;
	size_t pipeline_memory_limit;
	size_t user_requested_memory_for_allocation;
	bool is_init; /**< Flag to indicate whether cards on the machine have been checked/initialised. */
	aa_card_info selected_card_info; /** Stores all card information for all cards on the machine. */

	/** \brief Checks for GPUs on the machine.
	 * \details Provides the information back to the user via the card_info parameters.
	 * \returns A boolean to indicate whether this operation was successful.
	 */
	bool check_for_devices() {
		int deviceCount = 0;
		cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

		if (error_id != cudaSuccess) {
			LOG(log_level::dev_debug, "cudaGetDeviceCount returned" + std::to_string(error_id) + "->" + cudaGetErrorString(error_id));
			return false;
		}

		// This function call returns 0 if there are no CUDA capable devices.
		if (deviceCount == 0) {
			LOG(log_level::notice, "There are no available device(s) that support CUDA");
			return false;
		} else {
			LOG(log_level::notice, "Detected " + std::to_string(deviceCount) + " CUDA Capable device(s)");
		}

		if (deviceCount < selected_device_id) {
			LOG(log_level::notice, "Selected device with id " + std::to_string(selected_device_id) + " is not available");
			return false;
		}

		cudaSetDevice(selected_device_id);
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, selected_device_id);

		int driverVersion = 0;
		int runtimeVersion = 0;
		cudaDriverGetVersion(&driverVersion);
		cudaRuntimeGetVersion(&runtimeVersion);

		selected_card_info.name = deviceProp.name;
		selected_card_info.driver_version = driverVersion;
		selected_card_info.runtime_version = runtimeVersion;
		selected_card_info.global_memory = deviceProp.totalGlobalMem;
		selected_card_info.clock_rate = deviceProp.clockRate;
		selected_card_info.memory_clock_rate = deviceProp.memoryClockRate;
		selected_card_info.memory_bus_width = deviceProp.memoryBusWidth;
		selected_card_info.l2_cache_size = deviceProp.l2CacheSize;
		selected_card_info.compute_capability_major = deviceProp.major;
		selected_card_info.compute_capability_minor = deviceProp.minor;
		//All additional properties read out here...

		size_t free = 0;
		size_t total = 0;
		cudaMemGetInfo(&free, &total);
		selected_card_info.free_memory = free;
		selected_card_info.total_memory = total;
		LOG(log_level::notice, "Device info " + std::to_string(selected_device_id) + " free " + std::to_string(free) + " (" + selected_card_info.name.c_str() + ")");
		selected_card_info.user_requested_memory_for_allocation = 0;

		is_init = true;
		return true;
	}

	/** \returns A boolean to indicate whether selecting the card was successful. */
	bool init_card() {
		if (is_init==false) {
			return false;
		}

		#ifdef ASTRO_ACCELERATE_VERSION_H_DEFINED
		std::vector<int> compiled_cuda_sm_versions;
		std::stringstream s(ASTRO_ACCELERATE_CUDA_SM_VERSION);
		int i;
		while (s >> i) {
			compiled_cuda_sm_versions.push_back(i);

			if (s.peek() == ',') {
				s.ignore();
			}
		}
		#endif


		#ifdef ASTRO_ACCELERATE_VERSION_H_DEFINED
		bool found_compatible_architecture = false;
		std::string device_compute_capability = std::to_string(selected_card_info.compute_capability_major) + std::to_string(selected_card_info.compute_capability_minor);
		for (auto compiled_code : compiled_cuda_sm_versions) {
			if (std::to_string(compiled_code) > device_compute_capability) {
				continue;
			} else {
				LOG(log_level::notice, "Application binary compiled for compute capability "+std::to_string(compiled_code)+".");
				found_compatible_architecture = true;
			}
		}

		if (found_compatible_architecture) {
			LOG(log_level::notice, "The requested device has capability "+ device_compute_capability +" and the device code was compiled for a compatible architecture for the selected selected_card_info.");
		} else {
			LOG(log_level::error, "The requested device has capability "+ device_compute_capability +" but the device code was not compiled for a compatible architecture for the selected selected_card_info.");
			return false;
		}

		#else
		LOG(log_level::warning, "Because #include \"version.h\" is not created by this build system, the compute capability of the device code cannot be determined and compared with the physical device.");
		LOG(log_level::warning, "Please consider compiling using the CMakeLists file provided in the repository.");
		#endif

		return true;
	}

};

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DEVICE_INFO_HPP
