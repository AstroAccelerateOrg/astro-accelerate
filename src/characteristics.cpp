#include "presto.hpp"
#include "cdflib.hpp"

#include <stdio.h>
#include <stdlib.h>

namespace astroaccelerate {

static double extended_equiv_gaussian_sigma(double logp);
static double log_asymtotic_incomplete_gamma(double a, double z);
static double log_asymtotic_gamma(double z);




double extended_equiv_gaussian_sigma(double logp)
/*
  extended_equiv_gaussian_sigma(double logp):
      Return the equivalent gaussian sigma corresponding to the 
          natural log of the cumulative gaussian probability logp.
          In other words, return x, such that Q(x) = p, where Q(x)
          is the cumulative normal distribution.  This version uses
          the rational approximation from Abramowitz and Stegun,
          eqn 26.2.23.  Using the log(P) as input gives a much
          extended range.
*/
{
   double t, num, denom;

   t = sqrt(-2.0 * logp);
   num = 2.515517 + t * (0.802853 + t * 0.010328);
   denom = 1.0 + t * (1.432788 + t * (0.189269 + t * 0.001308));
   return t - num / denom;
}


double log_asymtotic_incomplete_gamma(double a, double z)
/*
  log_asymtotic_incomplete_gamma(double a, double z):
      Return the natural log of the incomplete gamma function in
          its asymtotic limit as z->infty.  This is from Abramowitz
          and Stegun eqn 6.5.32.
*/
{
   double x = 1.0, newxpart = 1.0, term = 1.0;
   int ii = 1;

   while (fabs(newxpart) > 1e-15) {
      term *= (a - ii);
      newxpart = term / pow(z, ii);
      x += newxpart;
      ii += 1;
   }
   return (a - 1.0) * log(z) - z + log(x);
}

double log_asymtotic_gamma(double z)
/*
  log_asymtotic_gamma(double z):
      Return the natural log of the gamma function in its asymtotic limit
          as z->infty.  This is from Abramowitz and Stegun eqn 6.1.41.
*/
{
   double x, y;

   x = (z - 0.5) * log(z) - z + 0.91893853320467267;
   y = 1.0 / (z * z);
   x += (((-5.9523809523809529e-4 * y
           + 7.9365079365079365079365e-4) * y
          - 2.7777777777777777777778e-3) * y + 8.3333333333333333333333e-2) / z;
   return x;
}


double equivalent_gaussian_sigma(double logp)
/* Return the approximate significance in Gaussian sigmas */
/* corresponding to a natural log probability logp        */
{
    double x;

    if (logp < -600.0) {
        x = extended_equiv_gaussian_sigma(logp);
    } else {
        int which, status;
        double p, q, bound, mean = 0.0, sd = 1.0;
        q = exp(logp);
        p = 1.0 - q;
        which = 2;
        status = 0;
        /* Convert to a sigma */
        cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound);
        if (status) {
            if (status == -2) {
                x = 0.0;
            } else if (status == -3) {
                x = 38.5;
            } else {
                printf("\nError in cdfnor() (candidate_sigma()):\n");
                printf("   status = %d, bound = %g\n", status, bound);
                printf("   p = %g, q = %g, x = %g, mean = %g, sd = %g\n\n",
                       p, q, x, mean, sd);
                exit(1);
            }
        }
    }
    if (x < 0.0)
        return 0.0;
    else
        return x;
}


double chi2_logp(double chi2, int dof)
/* Return the natural log probability corresponding to a chi^2 value */
/* of chi2 given dof degrees of freedom. */
{
    double logp;
    
    if (chi2 <= 0.0) {
        return -INFINITY;
    }

    if (chi2/dof > 15.0 || (dof > 150 && chi2/dof > 6.0)) {
      // printf("Using asymtotic expansion...\n");
        // Use some asymtotic expansions for the chi^2 distribution
        //   this is eqn 26.4.19 of A & S
        logp = log_asymtotic_incomplete_gamma(0.5*dof, 0.5*chi2) -
            log_asymtotic_gamma(0.5*dof);
    } else {
        int which, status;
        double p, q, bound, df = dof, x = chi2;

        which = 1;
        status = 0;
        /* Determine the basic probability */
        cdfchi(&which, &p, &q, &x, &df, &status, &bound);
        if (status) {
            printf("\nError in cdfchi() (chi2_logp()):\n");
            printf("   status = %d, bound = %g\n", status, bound);
            printf("   p = %g, q = %g, x = %g, df = %g\n\n",
                   p, q, x, df);
            exit(1);
        }
        // printf("p = %.3g  q = %.3g\n", p, q);
        logp = log(q);
    }
    return logp;
}


double chi2_sigma(double chi2, int dof)
/* Return the approximate significance in Gaussian sigmas        */
/* sigmas of a chi^2 value of chi2 given dof degrees of freedom. */
{
    double logp;

    if (chi2 <= 0.0) {
        return 0.0;
    }

    // Get the natural log probability
    logp = chi2_logp(chi2, dof);

    // Convert to sigma
    return equivalent_gaussian_sigma(logp);
}


double candidate_sigma(double power, int numsum, double numtrials)
/* Return the approximate significance in Gaussian       */
/* sigmas of a candidate of numsum summed powers,        */
/* taking into account the number of independent trials. */
{
    double logp, chi2, dof;

    if (power <= 0.0) {
        return 0.0;
    }

    // Get the natural log probability
    chi2 = 2.0 * power;
    dof = 2.0 * numsum;
    logp = chi2_logp(chi2, dof);

    // Correct for numtrials
    logp += log(numtrials);

    // Convert to sigma
    return equivalent_gaussian_sigma(logp);
}




/* Optional non-macro definitions of power and phase */
/*

__inline__ double POWER(double rl, double im)
{
  return rl * rl + im * im;
}


__inline__ double PHASE(double rl, double im)
{
  double temp;

  temp = 57.2957795130823208 * atan2(im, rl);
  return (temp > 0.0) ? temp : temp + 360.0;
}


double RADIAN_PHASE(double rl, double im)
{
  double temp;

  temp = atan2(im, rl);
  return (temp > 0.0) ? temp : temp + TWOPI;
}

*/

} //namespace astroaccelerate
