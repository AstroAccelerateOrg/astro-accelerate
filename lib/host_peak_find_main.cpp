#include <AstroAccelerate/device_peak_find.h>

#include "managed_allocator.hpp"

#include <fstream>
#include <iostream>
#include <vector>

#include <omp.h>

int main(int argc, char * argv[]) {
    std::cout << "Reading from input file: " << argv[1] << std::endl;
    std::ifstream infile;
    infile.open(argv[1], std::ios::in|std::ios::binary);
    infile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    if (infile.bad())
    {
        std::cerr << "Error opening file" << std::endl;
	exit(1);
    }

    infile.seekg(0, std::ios::end);
    auto sz = infile.tellg();
    infile.seekg(0, std::ios::beg);
    
    std::vector<float, managed_allocator<float>> data;
    data.reserve(sz/sizeof(float)/4);
    
    while (!infile.eof()) {
	float f1, f2, f3, f4;
	try {
            infile.read(reinterpret_cast<char*>(&f1), sizeof(float));
            infile.read(reinterpret_cast<char*>(&f2), sizeof(float));
            infile.read(reinterpret_cast<char*>(&f3), sizeof(float));
            infile.read(reinterpret_cast<char*>(&f4), sizeof(float));
        } catch (...) {
            break;
        }
	data.push_back(f3);
    }

    const int ntrials = 1;
    if (data.size()) {
	cudaSetDevice(0);
    	double start_t, end_t;
	auto width = data.size() / 4096;
        auto height = 4096;
        std::vector<unsigned short, managed_allocator<unsigned short>> output(data.size(), 0);
	peakfind(data.data(), 0, data.size(), 1, output.data());
	start_t = omp_get_wtime();
        for (int i=0; i < ntrials; ++i) {
	    peakfind(data.data(), 0, data.size(), 1, output.data());
	}
	cudaDeviceSynchronize();
	end_t = omp_get_wtime();
	float time = (float) ( end_t - start_t );
	auto data_mb = ntrials*data.size() * sizeof(float) / (1024 * 1024);
	std::cout << "Processed " <<  data_mb << " MB in " << time << " secs (" << data_mb/(1024*time) << ") GB/s" << std::endl;
	std::size_t peak_count = 0;
	for (int i=0; i < output.size() ; ++i) {
		if (output[i] == 1) {
		  ++peak_count;
		}
	}
	std::cout << "Peaks: " << peak_count <<  "/" << output.size() << std::endl;
        std::cout << "\n" << "Testing v1" << std::endl;
 	peakfind(data.data(), 0, width, height, output.data());
	start_t = omp_get_wtime();
        for (int i=0; i < ntrials; ++i) {
	    peakfind(data.data(), 0, width, height, output.data());
	}
	cudaDeviceSynchronize();
	end_t = omp_get_wtime();
	time = (float) ( end_t - start_t );
	data_mb = ntrials*data.size() * sizeof(float) / (1024 * 1024);
	std::cout << "Processed " <<  data_mb << " MB in " << time << " secs (" << data_mb/(1024*time) << ") GB/s" << std::endl;
	peak_count = 0;
	for (int i=0; i < output.size() ; ++i) {
		if (output[i] == 1) {
		  ++peak_count;
		}
	}
	std::cout << "Peaks: " << peak_count <<  "/" << output.size() << std::endl;
        std::cout << "\n" << "Testing v2" << std::endl;
	peakfind_v2(data.data(), 0, width, height, output.data());
	start_t = omp_get_wtime();
        for (int i=0; i < ntrials; ++i) {
	    peakfind_v2(data.data(), 0, width, height, output.data());
	}
	cudaDeviceSynchronize();
	end_t = omp_get_wtime();
	time = (float) ( end_t - start_t );
	data_mb = ntrials*data.size() * sizeof(float) / (1024 * 1024);
	std::cout << "Processed " <<  data_mb << " MB in " << time << " secs (" << data_mb/(1024*time) << ") GB/s" << std::endl;
	peak_count = 0;
	for (int i=0; i < output.size() ; ++i) {
		if (output[i] == 1) {
		  ++peak_count;
		  //std::cout << i << " " << data[i] << std::endl;
		}
	}
	std::cout << "Peaks: " << peak_count <<  "/" << output.size() << std::endl;
        std::cout << "\n" << "Testing v3" << std::endl;
	peakfind_v3(data.data(), 0, width, height, output.data());
	start_t = omp_get_wtime();
        for (int i=0; i < ntrials; ++i) {
	    peakfind_v3(data.data(), 0, width, height, output.data());
	}
	cudaDeviceSynchronize();
	end_t = omp_get_wtime();
	time = (float) ( end_t - start_t );
	data_mb = ntrials*data.size() * sizeof(float) / (1024 * 1024);
	std::cout << "Processed " <<  data_mb << " MB in " << time << " secs (" << data_mb/(1024*time) << ") GB/s" << std::endl;
	peak_count = 0;
	for (int i=0; i < output.size() ; ++i) {
		if (output[i] == 1) {
		  ++peak_count;
		  //std::cout << i << " " << data[i] << std::endl;
		}
	}
	std::cout << "Peaks: " << peak_count <<  "/" << output.size() << std::endl;
        std::cout << "\n" << "Testing v4" << std::endl;
	peakfind_v4(data.data(), 0, width, height, output.data());
	start_t = omp_get_wtime();
        for (int i=0; i < ntrials; ++i) {
	    peakfind_v4(data.data(), 0, width, height, output.data());
	}
	cudaDeviceSynchronize();
	end_t = omp_get_wtime();
	time = (float) ( end_t - start_t );
	data_mb = ntrials*data.size() * sizeof(float) / (1024 * 1024);
	std::cout << "Processed " <<  data_mb << " MB in " << time << " secs (" << data_mb/(1024*time) << ") GB/s" << std::endl;
	peak_count = 0;
	for (int i=0; i < output.size() ; ++i) {
		if (output[i] == 1) {
		  ++peak_count;
		  //std::cout << i << " " << data[i] << std::endl;
		}
	}
	std::cout << "Peaks: " << peak_count <<  "/" << output.size() << std::endl;
     } else {
        std::cerr << "No data to process!" << std::endl;
        exit(1);
    }
    return 0;
}
