#include <AstroAccelerate/device_peak_find.h>

#include "managed_allocator.hpp"

#include <fstream>
#include <iostream>
#include <vector>

#include <omp.h>

int main(int argc, char * argv[]) {
    std::vector<float, managed_allocator<float>> data;
    int ntrials = 10;
    if (argc > 1) {
        ntrials = atoi(argv[1]);
    }
    std::cerr << "Running with " << ntrials << " trials\n";

    while (true) {
	    data.resize(data.size() + 128ul*1024*1024/sizeof(float));
	    std::cout << "\n" << data.size() * sizeof(float) /1024 / 1024;
	    if (data.size()) {
		    cudaSetDevice(0);
		    double start_t, end_t;
		    auto width = data.size() / 4096;
		    auto height = 4096;
		    std::vector<unsigned short, managed_allocator<unsigned short>> output(data.size(), 0);
		    peakfind(data.data(), 0, width, height, output.data());
		    start_t = omp_get_wtime();
		    for (int i=0; i < ntrials; ++i) {
			    peakfind(data.data(), 0, width, height, output.data());
		    }
		    cudaDeviceSynchronize();
		    end_t = omp_get_wtime();
		    float time = (float) ( end_t - start_t );
		    auto data_mb = ntrials*data.size() * sizeof(float) / (1024 * 1024);
		    std::cout <<  " " << data_mb/(1024 * time) << std::flush;
		    std::size_t peak_count = 0;
		    for (int i=0; i < output.size() ; ++i) {
			    if (output[i] == 1) {
				    ++peak_count;
			    }
		    }
		    peakfind_v2(data.data(), 0, width, height, output.data());
		    start_t = omp_get_wtime();
		    for (int i=0; i < ntrials; ++i) {
			    peakfind_v2(data.data(), 0, width, height, output.data());
		    }
		    cudaDeviceSynchronize();
		    end_t = omp_get_wtime();
		    time = (float) ( end_t - start_t );
		    data_mb = ntrials*data.size() * sizeof(float) / (1024 * 1024);
		    std::cout <<  " " << data_mb/(1024 * time) << std::flush;
		    peak_count = 0;
		    for (int i=0; i < output.size() ; ++i) {
			    if (output[i] == 1) {
				    ++peak_count;
			    }
		    }
		    peakfind_v3(data.data(), 0, width, height, output.data());
		    start_t = omp_get_wtime();
		    for (int i=0; i < ntrials; ++i) {
			    peakfind_v3(data.data(), 0, width, height, output.data());
		    }
		    cudaDeviceSynchronize();
		    end_t = omp_get_wtime();
		    time = (float) ( end_t - start_t );
		    data_mb = ntrials*data.size() * sizeof(float) / (1024 * 1024);
		    std::cout <<  " " << data_mb/(1024 * time) << std::flush;
		    peak_count = 0;
		    for (int i=0; i < output.size() ; ++i) {
			    if (output[i] == 1) {
				    ++peak_count;
			    }
		    }
		    peakfind_v4(data.data(), 0, width, height, output.data());
		    start_t = omp_get_wtime();
		    for (int i=0; i < ntrials; ++i) {
			    peakfind_v4(data.data(), 0, width, height, output.data());
		    }
		    cudaDeviceSynchronize();
		    end_t = omp_get_wtime();
		    time = (float) ( end_t - start_t );
		    data_mb = ntrials*data.size() * sizeof(float) / (1024 * 1024);
		    std::cout <<  " " << data_mb/(1024 * time) << std::flush;
		    peak_count = 0;
		    for (int i=0; i < output.size() ; ++i) {
			    if (output[i] == 1) {
				    ++peak_count;
			    }
		    }
	    } else {
		    std::cerr << "No data to process!" << std::endl;
		    exit(1);
	    }
    }
    return 0;
}
