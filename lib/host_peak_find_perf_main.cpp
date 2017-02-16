#include <AstroAccelerate/device_peak_find.h>

#include "managed_allocator.hpp"
#include "accumulator.hpp"

#include <fstream>
#include <iostream>
#include <vector>

#include <omp.h>

template<typename Func>
void test_algorithm(const std::vector<float, managed_allocator<float>> & data, std::vector<unsigned short, managed_allocator<unsigned short>> output, int width, int height, Func peakfind, int ntrials, std::string algoname) {
    accumulator<float, double> stats;
    double start_t, end_t;
    for (int i=0; i < ntrials; ++i) {
	    start_t = omp_get_wtime();
	    peakfind(data.data(), 0, width, height, output.data(), 0.0f);
	    end_t = omp_get_wtime();
	    stats(end_t-start_t);
    }
    cudaDeviceSynchronize();
    std::size_t peak_count = 0;
    for (int i=0; i < output.size() ; ++i) {
	    if (output[i] == 1) {
		    ++peak_count;
	    }
    }
    auto time = stats.sum;
    auto data_mb = ntrials*data.size() * sizeof(float) / (1024 * 1024);
    std::cout << algoname << ", " << data.size() * sizeof(float) / (1024*1024) << ", " << data_mb/(1024 * time) << ", " << stats.mean() << ", " << stats.variance() << ", " << peak_count << std::endl;
}

int main(int argc, char * argv[]) {
    std::vector<float, managed_allocator<float>> data;
    int ntrials = 10;
    if (argc > 1) {
        ntrials = atoi(argv[1]);
    }
    std::cerr << "Running with " << ntrials << " trials\n";
    std::cout << "Algorithm, DataSize, data_rate(GB/s), mean time, variance, peaks\n";

    while (true) {
	    data.resize(data.size() + 128ul*1024*1024/sizeof(float));
	    if (data.size()) {
		    cudaSetDevice(0);
		    auto width = data.size() / 4096;
		    auto height = 4096;
		    std::vector<unsigned short, managed_allocator<unsigned short>> output(data.size(), 0);
		    peakfind(data.data(), 0, data.size(), 1, output.data());
                    test_algorithm(data, output, width, 1, peakfind, ntrials, "v1-1D");

		    peakfind(data.data(), 0, width, height, output.data());
                    test_algorithm(data, output, width, height, peakfind, ntrials, "v1");

		    peakfind_v2(data.data(), 0, width, height, output.data());
                    test_algorithm(data, output, width, height, peakfind_v2, ntrials, "v2");
		    
                    peakfind_v3(data.data(), 0, width, height, output.data());
                    test_algorithm(data, output, width, height, peakfind_v3, ntrials, "v3");
		    
                    peakfind_v4(data.data(), 0, width, height, output.data());
                    test_algorithm(data, output, width, height, peakfind_v4, ntrials, "v4");
	    } else {
		    std::cerr << "No data to process!" << std::endl;
		    exit(1);
	    }
    }
    return 0;
}
