#include <stdint.h>

/**
 * Finds the peaks in the 2d input array @p d_input.
 * @param d_input_linestep Defines the padding at the end of each row of the input data
 * @param width
 * @param height
 * @param d_output Output array - a value of 1 in this array indicates that it there is a peak in the corresponding input data. Any other values will be 0. Note that this memory is written to in a coalesced fashion so can be host mapped memory.
 */
void peakfind_npp(const float *d_input, int32_t d_input_linestep, int32_t width, int32_t height, unsigned short *d_output, const float threshold=0.0f);


/**
 * Finds the peaks in the 2d input array @p d_input.
 * @param d_input_linestep Defines the padding at the end of each row of the input data
 * @param width
 * @param height
 * @param d_output Output array - a value of 1 in this array indicates that it there is a peak in the corresponding input data. Any other values will be 0. Note that this memory is written to in a coalesced fashion so can be host mapped memory.
 */
void peakfind(const float *d_input, int32_t d_input_linestep, int32_t width, int32_t height, unsigned short *d_output, const float threshold=0.0f);

/**
 * Finds the peaks in the 2d input array @p d_input.
 * @param d_input_linestep Defines the padding at the end of each row of the input data
 * @param width
 * @param height
 * @param d_output Output array - a value of 1 in this array indicates that it there is a peak in the corresponding input data. Any other values will be 0. Note that this memory is written to in a coalesced fashion so can be host mapped memory.
 */
void peakfind_v2(const float *d_input, int32_t d_input_linestep, int32_t width, int32_t height, unsigned short *d_output, const float threshold=0.0f);

/**
 * Finds the peaks in the 2d input array @p d_input.
 * @param d_input_linestep Defines the padding at the end of each row of the input data
 * @param width
 * @param height
 * @param d_output Output array - a value of 1 in this array indicates that it there is a peak in the corresponding input data. Any other values will be 0. Note that this memory is written to in a coalesced fashion so can be host mapped memory.
 */
void peakfind_v3(const float *d_input, int32_t d_input_linestep, int32_t width, int32_t height, unsigned short *d_output, const float threshold=0.0f);

/**
 * Finds the peaks in the 2d input array @p d_input.
 * @param d_input_linestep Defines the padding at the end of each row of the input data
 * @param width
 * @param height
 * @param d_output Output array - a value of 1 in this array indicates that it there is a peak in the corresponding input data. Any other values will be 0. Note that this memory is written to in a coalesced fashion so can be host mapped memory.
 */
void peakfind_v4(const float *d_input, int32_t d_input_linestep, int32_t width, int32_t height, unsigned short *d_output, const float threshold=0.0f);
