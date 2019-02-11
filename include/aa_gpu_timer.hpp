#ifndef ASTRO_ACCELERATE_AA_GPU_TIMER_HPP
#define ASTRO_ACCELERATE_AA_GPU_TIMER_HPP

#include <cuda_runtime.h>

namespace astroaccelerate {

  /**
   * \struct aa_gpu_timer
   * \brief Wrapper around CUDA events which implements a timer functionality.
   */
  struct aa_gpu_timer {
    cudaEvent_t start;
    cudaEvent_t stop;

    aa_gpu_timer() {
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
    }

    ~aa_gpu_timer() {
      cudaEventDestroy(start);
      cudaEventDestroy(stop);
    }

    void Start() {
      cudaEventRecord(start, 0);
    }

    void Stop() {
      cudaEventRecord(stop, 0);
    }

    float Elapsed() {
      float elapsed = 0.0;
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&elapsed, start, stop);
      return elapsed;
    }
  };

} // namespace astroaccelerate

#endif  // ASTRO_ACCELERATE_AA_GPU_TIMER_HPP
