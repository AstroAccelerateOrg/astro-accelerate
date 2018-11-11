#ifndef ASTRO_ACCELERATE_PRESTO_FUNCS_HPP
#define ASTRO_ACCELERATE_PRESTO_FUNCS_HPP

#include "fresnl.hpp"
#include "median.hpp"
#include <cufft.h>

namespace astroaccelerate {

int presto_z_resp_halfwidth(double z, int accuracy);
cufftComplex *presto_gen_r_response(double roffset, int numbetween, int numkern);
cufftComplex *presto_gen_z_response( double z, int numkern, int numbetween);
void presto_place_complex_kernel(cufftComplex * kernel, int numkernel,
				 cufftComplex * result, int numresult);
void presto_dered_sig(cufftComplex * fft, int numamps);
void presto_norm(cufftComplex * fft, int numamps);

} //namespace astroaccelerate
  
#endif
