/************************** PRESTO accelsearch functions *************************/

#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include <math.h>
#include <string.h>

/* PRESTO defines */
#define NUMLOCPOWAVG  20 //Number of bins (total) to average for local power:  Must be an even number (1/2 on each side). 
#define DELTAAVGBINS  5 // Number of bins next to freq in question to ignore (on each side) when determining local power.
#define NUMFINTBINS   16 // Number of bins on each side of the central frequency to sum for Fourier interpolation (low accuracy)
#define HIGHACC 1 //Accuracy setting for calculating kernel half width
#define LOWACC 0

int fresnl(double xxa, double *ssa, double *cca);

float median(float arr[], int n);

int presto_z_resp_halfwidth(double z, int accuracy)
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a Fourier signal with constant f-dot. (i.e  */
  /*  a constant frequency derivative)                                 */
  /*  Arguments:                                                       */
  /*    'z' is the Fourier Frequency derivative (# of bins the signal  */
  /*       smears over during the observation).                        */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */
{
   int m;

   z = fabs(z);

   if (accuracy == HIGHACC) {
      m = (long) (z * (0.002057 * z + 0.0377) + NUMFINTBINS * 3);
      m += ((NUMLOCPOWAVG >> 1) + DELTAAVGBINS);

      /* Prevent the equation from blowing up in large z cases */

      if (z > 100 && m > 1.2 * z)
         m = 1.2 * z;

   } else {
      m = (long) (z * (0.00089 * z + 0.3131) + NUMFINTBINS);
      m = (m < NUMFINTBINS) ? NUMFINTBINS : m;

      /* Prevent the equation from blowing up in large z cases */

      if (z > 100 && m > 0.6 * z)
         m = 0.6 * z;
   }
   return m;
}



cufftComplex *presto_gen_r_response(double roffset, int numbetween, int numkern)
  /*  Generate a complex response function for Fourier interpolation.  */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 2 = interbins) */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */
{
   int ii;
   double tmp, sinc, s, c, alpha, beta, delta, startr, r;
   cufftComplex *response;

   /* Check that the arguments are OK */

   if (roffset < 0.0 || roffset >= 1.0) {
      printf("\n  roffset = %f (out of bounds) in gen_r_response().\n\n", roffset);
      exit(-1);
   }
   if (numbetween < 1 || numbetween >= 20000) {
      printf("\n  numbetween = %d (out of bounds) in gen_r_response().\n\n",
             numbetween);
      exit(-1);
   }
   if (numkern < numbetween) {
      printf("\n  numkern = %d (out of bounds) in gen_r_response().\n\n", numkern);
      exit(-1);
   }
   if ((numkern % (2 * numbetween)) != 0) {
      printf("\n  numkern %% (2 * numbetween) != 0 in gen_r_response().\n\n");
      exit(-1);
   }

   /* Prep the recursion */

   response = (cufftComplex*)malloc(numkern*sizeof(cufftComplex));          // gen_cvect(numkern);
   startr = M_PI * (numkern / (double) (2 * numbetween) + roffset);
   delta = -M_PI / numbetween;
   tmp = sin(0.5 * delta);
   alpha = -2.0 * tmp * tmp;
   beta = sin(delta);
   c = cos(startr);
   s = sin(startr);

   /* Generate the points */

   for (ii = 0, r = startr; ii < numkern; ii++, r += delta) {
      if (r == 0.0)
         sinc = 1.0;
      else
         sinc = s / r;
      response[ii].x = c * sinc;
      response[ii].y = s * sinc;
      c = alpha * (tmp = c) - beta * s + c;
      s = alpha * s + beta * tmp + s;
   }

   /* Correct for divide by zero when the roffset is close to zero */

   if (roffset < 1E-3) {
      response[numkern / 2].x = 1 - 6.579736267392905746 * (tmp = roffset * roffset);
      response[numkern / 2].y = roffset * (M_PI - 10.335425560099940058 * tmp);
   }
   return response;
}


cufftComplex *presto_gen_z_response( double z, int numkern, int numbetween)
  /*  Generate the response function for Fourier f-dot interpolation.  */
  /*  Arguments:                                                       */
  /*    'z' is the Fourier Frequency derivative (# of bins the signal  */
  /*       smears over during the observation).                        */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */
{
   int ii, signz, numkernby2;
   double absz, zd, tmp, r, xx, yy, zz, startr, startroffset;
   double fressy, frescy, fressz, frescz, tmprl, tmpim;
   double s, c, pibyz, cons, delta;
   cufftComplex *response;

   if (numbetween < 1 || numbetween >= 20000) {
      printf("\n  numbetween = %d (out of bounds) in gen_z_response().\n\n",
             numbetween);
      exit(-1);
   }
   if (numkern < numbetween) {
      printf("\n  numkern = %d (out of bounds) in gen_z_response().\n\n", numkern);
      exit(-1);
   }
   if ((numkern % (2 * numbetween)) != 0) {
      printf("\n  numkern %% (2 * numbetween) != 0 in gen_z_response().\n\n");
      exit(-1);
   }

   
   absz = fabs(z);
   if (absz < 1E-4) {
     response = presto_gen_r_response(0.0, numbetween, numkern);
     return response;
   }
   response = (cufftComplex*)malloc(numkern*sizeof(cufftComplex));          // gen_cvect(numkern);

   /* Begin the calculations */

   startr = - (0.5 * z);
   startroffset = (startr < 0) ? 1.0 + modf(startr, &tmprl) : modf(startr, &tmprl);
   signz = (z < 0.0) ? -1 : 1;
   zd = signz * M_SQRT2 / sqrt(absz);
   cons = zd / 2.0;
   pibyz = M_PI / z;
   startr += numkern / (2.0*numbetween);
   delta = -1.0 / numbetween ;

   for (ii = 0, r = startr; ii < numkern; ii++, r += delta) {
      yy = r * zd;
      zz = yy + z * zd;
      xx = pibyz * r * r;
      c = cos(xx);
      s = sin(xx);
      fresnl(yy, &fressy, &frescy);
      fresnl(zz, &fressz, &frescz);
      tmprl = signz * (frescz - frescy);
      tmpim = fressy - fressz;
      response[ii].x = ((tmp = tmprl) * c - tmpim * s) * cons;
      response[ii].y = -(tmp * s + tmpim * c) * cons;
   }

   /* Correct for divide by zero when the roffset and z is close to zero */

   if (startroffset < 1E-3 && absz < 1E-3) {
      zz = z * z;
      xx = startroffset * startroffset;
      numkernby2 = numkern / 2;
      response[numkernby2].x = 1.0 - 0.16449340668482264365 * zz;
      response[numkernby2].y = -0.5235987755982988731 * z;
      response[numkernby2].x += startroffset * 1.6449340668482264365 * z;
      response[numkernby2].y += startroffset * (M_PI - 0.5167712780049970029 * zz);
      response[numkernby2].x += xx * (-6.579736267392905746
                                      + 0.9277056288952613070 * zz);
      response[numkernby2].y += xx * (3.1006276680299820175 * z);
   }
   int j;
   for ( j=0; j<numkern; j++){
     if(isnan(response[j].x)||isnan(response[j].y)){
       printf("\nnan detected in template in point %d: %f\t%f\n", j, response[j].x, response[j].y);
       //      exit(1);
     }
   }


   return response;
}


void presto_place_complex_kernel(cufftComplex * kernel, int numkernel,
                          cufftComplex * result, int numresult)
  /* This routine places the kernel in a zero filled array */
  /* with half of the response at the beginning and half   */
  /* of the response at the end of the result array.  See  */
  /* Numerical Recipes in C 2ed, p 541 for more info.      */
  /* Arguments:                                            */
  /*   'kernel' is a complex response function.  Bin zero  */
  /*      response is in bin numkernel/2.                  */
  /*   'numkernel' is the number of points in the kernel.  */
  /*      This should be an even number.                   */
  /*   'result' is the result array.                       */
  /*   'numresult' is the number of points in the result.  */
{
  int ii, halfwidth;
  cufftComplex zeros;
  zeros.x = zeros.y = 0.0;
  
  halfwidth = numkernel / 2;
  
  for (ii = 0; ii < numresult; ii++){
    //    printf("\nstart\n");
    result[ii] = zeros;
  }
  //printf("\nmemcpy numkernel = %d  halfwidth=%d numresult=%d\n",numkernel,halfwidth,numresult);
  memcpy(result, kernel + halfwidth, sizeof(cufftComplex) * halfwidth);
  //  printf("\nmemcpy 1st half ok\n");
  memcpy(result + numresult - halfwidth, kernel, sizeof(cufftComplex) * halfwidth);
  //  printf("\nmemcpy 2nd half ok\n");
}

void presto_dered_sig(cufftComplex * fft, int numamps)
/* Attempt to remove rednoise from a time series by using   */
/* a median-filter of logarithmically increasing width.     */
/* Thanks to Jason Hessels and Maggie Livingstone for the   */
/* initial implementation (in rednoise.c)                   */
{
  int ii, initialbuflen = 6, lastbuflen, newbuflen, maxbuflen =200;
  int newoffset = 1, fixedoffset = 1;
  float *powers, mean_old, mean_new;
  double slope, lineval, scaleval = 1.0, lineoffset = 0;

  powers = (float*)malloc(numamps*sizeof(float));

  /* Takes care of the DC term */
  fft[0].x = 1.0f;
  fft[0].y = 0.0f;

  /* Step through the input FFT and create powers */
  for (ii = 0; ii < numamps; ii++) {
    float powargr, powargi;
    powers[ii] = fft[ii].x*fft[ii].x + fft[ii].y*fft[ii].y;
  }

  /* Calculate initial values */
  mean_old = median(powers + newoffset, initialbuflen) / log(2.0);
  newoffset += initialbuflen;
  lastbuflen = initialbuflen;
  newbuflen = initialbuflen * log(newoffset);
  if (newbuflen > maxbuflen)
    newbuflen = maxbuflen;

  while (newoffset + newbuflen < numamps) {
    /* Calculate the next mean */
    mean_new = median(powers + newoffset, newbuflen) / log(2.0);
    //slope = (mean_new-mean_old)/(0.5*(newbuflen+lastbuflen));
    slope = (mean_new - mean_old) / (newbuflen + lastbuflen);

    /* Correct the previous segment */
    lineoffset = 0.5 * (newbuflen + lastbuflen);
    for (ii = 0; ii < lastbuflen; ii++) {
      lineval = mean_old + slope * (lineoffset - ii);
      scaleval = 1.0 / sqrt(lineval);
      fft[fixedoffset + ii].x *= scaleval;
      fft[fixedoffset + ii].y *= scaleval;
    }

    /* Update our values */
    fixedoffset += lastbuflen;
    lastbuflen = newbuflen;
    mean_old = mean_new;
    newoffset += lastbuflen;
    newbuflen = initialbuflen * log(newoffset);
    if (newbuflen > maxbuflen)
      newbuflen = maxbuflen;
  }

  /* Scale the last (partial) chunk the same way as the last point */

  while (fixedoffset < numamps) {
    fft[fixedoffset].x *= scaleval;
    fft[fixedoffset].y *= scaleval;
    fixedoffset++;
  }

  /* Free the powers */
  free(powers);
}


void presto_norm(cufftComplex * fft, int numamps)
{
  float *powers;
  double norm;
  int ii;
  powers = (float*)malloc(numamps*sizeof(float));
  //fft[0].x=1.0f;
  //fft[0].y=1.0f;
  /* Step through the input FFT and create powers */
  for (ii = 0; ii < numamps; ii++) {
    powers[ii] = fft[ii].x*fft[ii].x + fft[ii].y*fft[ii].y;
  }

  norm = 1.0 / sqrt(median(powers, numamps)/log(2.0));
  
  for (ii = 0; ii < numamps; ii++) {
    fft[ii].x *= norm;
    fft[ii].y *= norm;
  }

free(powers);
}
