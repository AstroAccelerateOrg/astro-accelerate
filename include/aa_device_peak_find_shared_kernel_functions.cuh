#ifndef ASTRO_ACCELERATE_AA_DEVICE_PEAK_FIND_SHARED_KERNEL_FUNCTIONS_CUH
#define ASTRO_ACCELERATE_AA_DEVICE_PEAK_FIND_SHARED_KERNEL_FUNCTIONS_CUH

namespace astroaccelerate {

struct float3x3 {
  float x1, x2, x3;
  float y1, y2, y3;
  float z1, z2, z3;
};

__device__ __inline__ ushort is_peak(const float original_value, const float dilated_value, const float threshold) {
  return (original_value > threshold && original_value == dilated_value) ? 1u: 0u;
}

__device__ __inline__ float3 load_row(const float *input, int idx) {
  float3 val;
  val.x = __ldg(input+idx);
  val.y = __ldg(input+idx+1);
  val.z = __ldg(input+idx+2);
  return val;
}

__device__ __inline__ float2 load_row2(const float *input, int idx) {
  float2 val;
  val.x=__ldg(input+idx);
  val.y=__ldg(input+idx+1);
  return val;
}

__device__ __inline__ float3x3 load_block(const float *input, int idxX, int idxY, int width) {
  int idx = idxY*width + idxX;
  float3x3 data;
  float3 yrow = load_row(input, idx-1);
  float3 xrow = load_row(input, idx-width-1);
  float3 zrow = load_row(input, idx+width-1);
  data.x1 = xrow.x;
  data.x2 = xrow.y;
  data.x3 = xrow.z;
  data.y1 = yrow.x;
  data.y2 = yrow.y;
  data.y3 = yrow.z;
  data.z1 = zrow.x;
  data.z2 = zrow.y;
  data.z3 = zrow.z;
  return data;
}

__device__ __inline__ float3x3 load_block_top(const float *input, int idxX, int idxY, int width) {
  int idx = idxY*width + idxX;
  float3x3 data;
  float3 yrow = load_row(input, idx-1);
  float3 zrow = load_row(input, idx+width-1);
  data.y1 = yrow.x;
  data.y2 = yrow.y;
  data.y3 = yrow.z;
  data.z1 = zrow.x;
  data.z2 = zrow.y;
  data.z3 = zrow.z;
  return data;
}

__device__ __inline__ float3x3 load_block_bottom(const float *input, int idxX, int idxY, int width) {
  int idx = idxY*width + idxX;
  float3x3 data;
  float3 yrow = load_row(input, idx-1);
  float3 xrow = load_row(input, idx-width-1);
  data.y1 = yrow.x;
  data.y2 = yrow.y;
  data.y3 = yrow.z;
  data.x1 = xrow.x;
  data.x2 = xrow.y;
  data.x3 = xrow.z;
  return data;
}

__device__ __inline__ float3x3 load_block_left(const float *input, int idxX, int idxY, int width) {
  int idx = idxY*width + idxX;
  float3x3 data;
  float2 xrow = load_row2(input, idx-width);
  float2 yrow = load_row2(input, idx);
  float2 zrow = load_row2(input, idx+width);
  data.x2 = xrow.x;
  data.x3 = xrow.y;
  data.y2 = yrow.x;
  data.y3 = yrow.y;
  data.z2 = zrow.x;
  data.z3 = zrow.y;
  return data;
}

__device__ __inline__ float3x3 load_block_right(const float *input, int idxX, int idxY, int width) {
  int idx = idxY*width + idxX;
  float3x3 data;
  float2 xrow = load_row2(input, idx-width-1);
  float2 yrow = load_row2(input, idx-1);
  float2 zrow = load_row2(input, idx+width-1);
  data.x1 = xrow.x;
  data.x2 = xrow.y;
  data.y1 = yrow.x;
  data.y2 = yrow.y;
  data.z1 = zrow.x;
  data.z2 = zrow.y;
  return data;
}

__device__ __inline__ float4 load_block_2x2(const float *input, int width) {
  float2 first, second;
  first.x  = __ldg(input);        first.y  = __ldg(input+1);
  second.x = __ldg(input+width);  second.y = __ldg(input+width+1);
  float4 result;
  result.x = first.x; result.y = first.y; result.z = second.x; result.w = second.y;
  return result;
}

__device__ __inline__ float dilate3x3_left(const float3x3 i) {
  float max = fmaxf(i.x2, i.y2);
  max = fmaxf(max, i.x3);
  max = fmaxf(max, i.y3);
  max = fmaxf(max, i.z2);
  max = fmaxf(max, i.z3);
  return max;
}

__device__ __inline__ float dilate3x3_right(const float3x3 i) {
  float max = fmaxf(i.x2, i.y2);
  max = fmaxf(max, i.x1);
  max = fmaxf(max, i.y1);
  max = fmaxf(max, i.z2);
  max = fmaxf(max, i.z1);
  return max;
}

__device__ __inline__ float dilate3x3_top(const float3x3 i) {
  float max = fmaxf(i.y1, i.y2);
  max = fmaxf(max, i.y3);
  max = fmaxf(max, i.z1);
  max = fmaxf(max, i.z2);
  max = fmaxf(max, i.z3);
  return max;
}

__device__ __inline__ float dilate3x3_bottom(const float3x3 i) {
  float max = fmaxf(i.y1, i.y2);
  max = fmaxf(max, i.y3);
  max = fmaxf(max, i.x1);
  max = fmaxf(max, i.x2);
  max = fmaxf(max, i.x3);
  return max;
}

__device__ __inline__ float dilate4(const float4 i) {
  float max = fmaxf(i.x, i.y);
  max = fmaxf(max, i.z);
  max = fmaxf(max, i.w);
  return max;
}

__device__ __inline__ float dilate3x3(const float3x3 i) {
  float max = fmaxf(i.x1, i.x2);
  max = fmaxf(max, i.x3);
  max = fmaxf(max, i.y1);
  max = fmaxf(max, i.y2);
  max = fmaxf(max, i.y3);
  max = fmaxf(max, i.z1);
  max = fmaxf(max, i.z2);
  max = fmaxf(max, i.z3);
  return max;
}

} // namespace astroaccelerate
#endif // ASTRO_ACCELERATE_AA_DEVICE_PEAK_FIND_SHARED_KERNEL_FUNCTIONS_CUH
