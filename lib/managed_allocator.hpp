#pragma once

#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>

template<class T>
class managed_allocator
{
  public:
    using value_type = T;
  
    value_type* allocate(size_t n)
    {
      value_type* result = nullptr;
  
      cudaError_t error = cudaMallocManaged(reinterpret_cast<void**>(&result), n*sizeof(T), cudaMemAttachGlobal);
  
      if(error != cudaSuccess)
      {
        throw thrust::system_error(error, thrust::cuda_category(), "managed_allocator::allocate(): cudaMallocManaged");
      }
  
      return result;
    }
  
    void deallocate(value_type* ptr, size_t)
    {
      cudaError_t error = cudaFree(ptr);
  
      if(error != cudaSuccess)
      {
        throw thrust::system_error(error, thrust::cuda_category(), "managed_allocator::deallocate(): cudaFree");
      }
    }
};

