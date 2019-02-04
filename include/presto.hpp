/* PRESTO Function declarations */
#ifndef ASTRO_ACCELERATE_PRESTO_HPP
#define ASTRO_ACCELERATE_PRESTO_HPP

#include <math.h>
#include <cufft.h>

namespace astroaccelerate {

#define LOWACC 0
#define SLIGHT 299792458.0
#ifndef SQRT2
#define SQRT2         1.4142135623730950488016887242096980785696718753769
#endif
#ifndef DBLCORRECT
#define DBLCORRECT    1e-14
#endif
#ifndef PI
#define PI            3.1415926535897932384626433832795028841971693993751
#endif
#ifndef TWOPI
#define TWOPI         6.2831853071795864769252867665590057683943387987502
#endif
#ifndef DEGTORAD
#define DEGTORAD      0.017453292519943295769236907684886127134428718885417
#endif
#ifndef RADTODEG
#define RADTODEG      57.29577951308232087679815481410517033240547246656
#endif
#ifndef PIBYTWO
#define PIBYTWO       1.5707963267948966192313216916397514420985846996876
#endif
#ifndef SOL
#define SOL           299792458.0
#endif
#ifndef SECPERJULYR
#define SECPERJULYR   31557600.0
#endif
#ifndef SECPERDAY
#define SECPERDAY     86400.0
#endif
#ifndef ARCSEC2RAD
#define ARCSEC2RAD    4.8481368110953599358991410235794797595635330237270e-6
#endif
#ifndef SEC2RAD
#define SEC2RAD       7.2722052166430399038487115353692196393452995355905e-5
#endif
#ifndef __GNUC__
#define __inline__
#endif

/* Maximum number of input files to try and patch together */
#define MAXPATCHFILES 100
/* Blocksize to use when reading datafiles or subbands */
#define SUBSBLOCKLEN 1024

/* various function-like macros */

#ifndef SWAP
/* Swaps two variables of undetermined type */
#define SWAP(a,b) tempzz=(a);(a)=(b);(b)=tempzz;
#endif

#ifndef POWER
/* Returns unnormalized Fourier power  */
/*   Requires the following variables in calling function */
/*   double powargr, powargi; */
#define POWER(r,i) (powargr=(r),powargi=(i),\
		    powargr*powargr+powargi*powargi)
#endif

#ifndef PHASE
/*   Returns Fourier phase (degrees)  */
/*   Requires the following variables in calling function */
/*   double phsargr, phsargi, phstmp; */
#define PHASE(r,i) (phsargr=(r),phsargi=(i),\
		    ((phstmp=RADTODEG*atan2(phsargi,phsargr)) > 0.0) ? \
		    phstmp : phstmp+360.0)
#endif
  
#ifndef RADIAN_PHASE
/* Returns Fourier phase (radians)  */
/*   Requires the following variables in calling function */
/*   double radargr, radargi, radtmp; */
#define RADIAN_PHASE(r,i) (radargr=(r),radargi=(i),\
		    ((radtmp=atan2(radargi,radargr)) > 0.0) ? \
		    radtmp : radtmp+TWOPI)
#endif

#define GET_BIT(c, n) (*(c+(n>>3)) >> (7-(n&7)) & 1)
#define SET_BIT(c, n) (*(c+(n>>3)) |= 1 << (7-(n&7)))
#define UNSET_BIT(c, n) (*(c+(n>>3)) &= ~(1 << (7-(n&7))))


//extern "C" int fresnl(double xxa, double *ssa, double *cca);

 //extern "C" 
cufftComplex *gen_z_response( double z, int numkern, int numbetween);

//extern "C" 
void place_complex_kernel(cufftComplex * kernel, int numkernel, cufftComplex * result, int numresult);

//extern "C" 
int z_resp_halfwidth(double z, int accuracy);

//extern "C" float median(float arr[], int n);

//extern "C" 
void dereddensig(cufftComplex * fft, int numamps);

//extern "C" 
void presto_norm(cufftComplex * fft, int numamps);

//extern "C" 
double equivalent_gaussian_sigma(double logp);
/* Return the approximate significance in Gaussian sigmas */
/* corresponding to a natural log probability logp        */

//extern "C" 
double chi2_logp(double chi2, int dof);
/* Return the natural log probability corresponding to a chi^2 value */
/* of chi2 given dof degrees of freedom. */

//extern "C" 
double chi2_sigma(double chi2, int dof);
/* Return the approximate significance in Gaussian sigmas        */
/* sigmas of a chi^2 value of chi2 given dof degrees of freedom. */

//extern "C" 
double candidate_sigma(double power, int numsum, double numtrials);
/* Return the approximate significance in Gaussian       */
/* sigmas of a candidate of numsum summed powers,        */
/* taking into account the number of independent trials. */

/* Return the natural log probability corresponding to a chi^2 value */
/* of chi2 given dof degrees of freedom. */

} // namespace astroaccelerate
  
#endif /* PRESTO_H */
