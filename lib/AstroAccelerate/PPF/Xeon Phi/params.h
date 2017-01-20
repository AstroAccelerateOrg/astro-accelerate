// RpCl is registers per Cache-line. It says how many registers one need to load a cache line. In case of Xeon Phi it is 1.
#define RpCl 1
// FpR is number of single precision numbers per AVX register. In case of Xeon Phi it is 512bit/32bit=16
#define FpR 16
// CpB is Cache-lines per block, where block is meant the division of data block. Cache line is 64bytes which means two AVX registers (256bit/32bit)=8, so one AVX register can hold 8 single precision numbers. Cache line is 64/4=16 single precision numbers. Thus CpB gives one dimension of the block C_B in our article.
#define CpB 16
// SpT is Spectra per thread, it says how many filtered spectra are produced by one thread. It defines second dimension of our grid (S_B in article)
#define SpT 7
// NUMTHREADS says how many threads code will use. The thread affinity has also effect on result and number of threads is only half of the story, sort of...
#define NUMTHREADS 240

// Note: When we have performed our auto-tuning we also changed value RpCl. Although this is bit dirty, since increasing RpCl from 1 to 2 essentially means multiplying CpB by two.
// For example RpCl 4
