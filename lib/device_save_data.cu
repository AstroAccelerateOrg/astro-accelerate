//{{{ save_data_from_device_to_host

//void save_data(cudaStream_t stream, float *device_pointer, float *host_pointer, size_t size) {
void save_data(float *device_pointer, float *host_pointer, size_t size)
{

	//{{{ Copy data and set up the GPU constants/variables.

//	cudaMemcpy(host_pointer, device_pointer, size, cudaMemcpyDeviceToHost);
	//cudaMemcpyAsync(host_pointer, device_pointer, size, cudaMemcpyDeviceToHost, stream);

	//}}}

}

//}}}
