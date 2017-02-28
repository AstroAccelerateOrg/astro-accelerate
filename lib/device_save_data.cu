//{{{ save_data_from_device_to_host

//void save_data(cudaStream_t stream, float *device_pointer, float *host_pointer, size_t size) {
void save_data(float *device_pointer, float *host_pointer, size_t size)
{

	//{{{ Copy data and set up the GPU constants/variables.

//	cudaMemcpy(host_pointer, device_pointer, size, cudaMemcpyDeviceToHost);
	//cudaMemcpyAsync(host_pointer, device_pointer, size, cudaMemcpyDeviceToHost, stream);

	//}}}

}

//void save_data(cudaStream_t stream, float *device_pointer, float *host_pointer, size_t size) {
void save_data_offset(float *device_pointer, int device_offset, float *host_pointer, int host_offset, size_t size)
{

	//{{{ Copy data and set up the GPU constants/variables.

	cudaMemcpy(host_pointer + host_offset, device_pointer + device_offset, size, cudaMemcpyDeviceToHost);
	//cudaMemcpyAsync(host_pointer, device_pointer, size, cudaMemcpyDeviceToHost, stream);

	//}}}

}
//}}}
