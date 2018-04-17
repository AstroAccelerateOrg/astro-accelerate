#ifndef GPU_TIMER_H__
#define GPU_TIMER_H__

#include <cuda_runtime.h>

struct GpuTimer
{
  cudaEvent_t start;
  cudaEvent_t stop;

  GpuTimer()
  {
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
  }

  ~GpuTimer()
  {
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
  }

  void Start()
  {
    cudaEventRecord(start, 0);
  }

  void StartWithStream(cudaStream_t stream)
  {
    cudaEventRecord(start, stream);
  }

  void Stop()
  {
    cudaEventRecord(stop, 0);
  }

  void StopWithStream(cudaStream_t stream)
  {
    cudaEventRecord(stop, stream);
  }

  float Elapsed()
  {
    float elapsed;
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed, start, stop);
    return elapsed;
  }

  float ElapsedWithStream(cudaStream_t stream)
	  {
    float elapsed;
	bool process_stop = false;
	while (process_stop == false) 
		{
			process_stop = cudaEventQuery(stop) == cudaSuccess;
		}
//    cudaStreamWaitEvent(stream, stop, 0);
//    cudaStreamWaitEvent(stream, stop, 0);
    cudaEventElapsedTime(&elapsed, start, stop);
    return elapsed;
  }

};

#endif  /* GPU_TIMER_H__ */
