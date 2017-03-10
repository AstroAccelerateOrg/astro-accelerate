// Warning: Calculation of these parameters is different for each precision case. Sorry.

// This is dummy parameter, which is not used in the code itself. It says how many thread-blocks we would like to be resident on single SM (multiprocessor)
#define ACTIVE_BLOCKS 3
// This is again dummy parameter. It says how many single precision floating point numbers can fit into shared memory available to single block
#define TOT_SM_SIZE 12288
// This is again dummy parameter which says how many channels are processed per single thread-block. It is accualy size of a warp=32.
#define CHANNELS_PER_BLOCK 32
// note: for Maxwell generation the ACTIVE_BLOCKS could be half of blocks present on single SM as this generation has 96kB of shared memory, but shared memory per block is still 48kB

// gives number of taps
#define TAPS 16
// COEFF_SIZE is given by number of taps and channels processed per thread-block COEFF_SIZE=CHANNELS_PER_BLOCK*TAPS=512
#define COEFF_SIZE 512
// DATA_SIZE says how many input data elements (in floating point numbers) we want to store in shared memory per thread-block. DATA_SIZE=TOT_SM_SIZE/ACTIVE_BLOCKS=4096; DATA_SIZE=DATA_SIZE-COEFF_SIZE=3584; this is because we need to store coefficients in the shared memory along with the input data. Lastly we must divide this by two DATA_SIZE=DATA_SIZE/2=1792; This is because we store real and imaginary input data separately to prevent bank conflicts. 
#define DATA_SIZE 1792
// THREADS_PER_BLOCK gives number of threads per thread-block. It could be calculated as such THREADS_PER_BLOCK=MAX_THREADS_PER_SM/ACTIVE_BLOCKS; rounded to nearest lower multiple of 32; In case of Maxwell generation MAX_THREADS_PER_SM=2048, thus THREADS_PER_BLOCK=682.6666, rounding to nearest lower multiple of 32 gives THREADS_PER_BLOCK=672;
#define THREADS_PER_BLOCK 672
// SUBBLOCK_SIZE gives size of the sub-block as given in our article. It is calculated as ratio of DATA_SIZE and THREADS_PER_BLOCK rounded up so SUBBLOCK_SIZE=DATA_SIZE/THREADS_PER_BLOCK=2.6666, rounding up gives SUBBLOCK_SIZE=3;
#define SUBBLOCK_SIZE 3