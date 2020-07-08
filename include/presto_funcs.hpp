#ifndef ASTRO_ACCELERATE_PRESTO_FUNCS_HPP
#define ASTRO_ACCELERATE_PRESTO_FUNCS_HPP

#include <cufft.h>

#include "fresnl.hpp"
#include "aa_median.hpp"

#include "aa_jerk_plan.hpp"
//#include "aa_jerk_strategy.hpp"

namespace astroaccelerate {

	int presto_r_resp_halfwidth(int accuracy);
	int presto_z_resp_halfwidth(double z, int accuracy);
	int presto_w_resp_halfwidth(double z, double w, int accuracy);
	

	cufftComplex *presto_gen_r_response(double roffset, int numbetween, int numkern);
	cufftComplex* presto_gen_z_response(double roffset, int numbetween, double z, int numkern);
	cufftComplex *presto_gen_w_response(double roffset, int numbetween, double z, double w, int numkern);
	
	void presto_place_complex_kernel(cufftComplex * kernel, int numkernel, cufftComplex * result, int numresult);
	void presto_dered_sig(cufftComplex * fft, int numamps);
	void presto_norm(cufftComplex * fft, int numamps);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_PRESTO_FUNCS_HPP
