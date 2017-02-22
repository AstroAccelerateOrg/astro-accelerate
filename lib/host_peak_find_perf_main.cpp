#include <AstroAccelerate/device_peak_find.h>

#include "managed_allocator.hpp"
#include "accumulator.hpp"

#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <condition_variable>

#include <unistd.h>
#include <omp.h>

#include "nvml.h"

class PowerMonitor {
public:
    PowerMonitor() : active(true), pm(&PowerMonitor::monitor_power, this), max_power(0) {
        sleep(1);
    }

    ~PowerMonitor() {
        active = false;
        m_end = true;
        start_cond.notify_all();

        pm.join();
    }

    void wait_for_start() {
        std::unique_lock<std::mutex> lock(m);
        m_end = false;
        start_cond.wait(lock);
    }

    void start() {
        start_cond.notify_all();
    }

    void end() {
        m_end = true; 
    }

    void csv_report() {
        std::cout << stats.N << ", " << max_power << ", " << stats.mean() << ", " << stats.variance();
    }

private:
    void monitor_power() {
         nvmlDevice_t device;
         unsigned int power;
         nvmlDeviceGetHandleByIndex(0, &device);
         while (active) {
             wait_for_start();
             m.unlock();
             do {
                 nvmlDeviceGetPowerUsage(device, &power);
                 max_power = std::max(max_power, power);
                 stats(power);
                 usleep(1);
             } while (!m_end);
         }
    }

private:
    volatile bool active; //Has to be above the thread so its initialised first
    std::thread pm;

    std::mutex m;
    std::condition_variable start_cond;
    volatile bool m_end;
    accumulator<unsigned int, unsigned long long> stats;
    unsigned int max_power;
};

template<typename Func>
void test_algorithm(const std::vector<float, managed_allocator<float>> & data, std::vector<unsigned short, managed_allocator<unsigned short>> output, int width, int height, Func peakfind, int ntrials, std::string algoname) {
    PowerMonitor power_monitor;
    accumulator<float, double> stats;
    double start_t, end_t;
    for (int i=0; i < ntrials; ++i) {
            power_monitor.start();
	    start_t = omp_get_wtime();
	    peakfind(data.data(), 0, width, height, output.data(), 0.0f);
	    end_t = omp_get_wtime();
            power_monitor.end();
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
    auto data_gb = ntrials*data.size() * sizeof(float) / (1024 * 1024 * 1024);
    std::cout << algoname << ", " << data.size() * sizeof(float) / (1024*1024) << ", " << data_gb/time << ", " << stats.mean() << ", " << stats.variance() << ", " << peak_count << ", ";
    power_monitor.csv_report();
    std::cout << std::endl;
}

int main(int argc, char * argv[]) {
    nvmlInit();
    std::vector<float, managed_allocator<float>> data;
    int ntrials = 10;
    if (argc > 1) {
        ntrials = atoi(argv[1]);
    }
    std::cerr << "Running with " << ntrials << " trials\n";
    std::cout << "Algorithm, DataSize(MB), data_rate(GB/s), mean time, variance, peaks, power samples, max power (mW), mean power (mW), power variance\n";

    while (true) {
	    data.resize(data.size() + 128ul*1024*1024/sizeof(float));
	    if (data.size()) {
		    cudaSetDevice(0);
		    auto width = data.size() / 4096;
		    auto height = 4096;
		    std::vector<unsigned short, managed_allocator<unsigned short>> output(data.size(), 0);
		    peakfind(data.data(), 0, data.size(), 1, output.data());
                    test_algorithm(data, output, data.size(), 1, peakfind, ntrials, "v1-1D");

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
