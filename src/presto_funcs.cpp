/************************** PRESTO accelsearch functions *************************/

#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include <string.h>
#include "fresnl.hpp"
#include "presto_funcs.hpp"
#include "aa_jerk_plan.hpp"
#include "aa_jerk_strategy.hpp"

namespace astroaccelerate {

  /* PRESTO defines */
#define NUMLOCPOWAVG  20 //Number of bins (total) to average for local power:  Must be an even number (1/2 on each side). 
#define DELTAAVGBINS  5 // Number of bins next to freq in question to ignore (on each side) when determining local power.
#define NUMFINTBINS   16 // Number of bins on each side of the central frequency to sum for Fourier interpolation (low accuracy)
#define HIGHACC 1 //Accuracy setting for calculating kernel half width
#define LOWACC 0
#define TWOPI 6.283185307179586476925286766559

  int fresnl(double xxa, double *ssa, double *cca);

  float median(float arr[], int n);

long long next2_to_n(long long x)
/* Return the first value of 2^n >= x */
{
    long long i = 1;

    while (i < x)
        i <<= 1;
    return i;
}

int presto_r_resp_halfwidth(int accuracy){
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a standard Fourier signal.                  */
  /*  Arguments:                                                       */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */
    if (accuracy == HIGHACC) {
        return ((NUMFINTBINS * 3) + (NUMLOCPOWAVG >> 1) + DELTAAVGBINS);
    } else {
        return NUMFINTBINS;
    }
}


int presto_z_resp_halfwidth(double z, int accuracy) {
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
    int m = (int) (0.5 * 1.1 * fabs(z));
    if (accuracy == HIGHACC) {
        m += NUMFINTBINS * 3;
    } else {
        m += NUMFINTBINS;
    }
    return m;
}


int presto_w_resp_halfwidth(double z, double w, int accuracy){
  /*  Return the approximate kernel half width in FFT bins required    */
  /*  to achieve a fairly high accuracy correlation based correction   */
  /*  or interpolation for a Fourier signal with an f-dot that (i.e    */
  /*  varies linearly in time -- a constant f-dotdot)                  */
  /*  Arguments:                                                       */
  /*    'z' is the average Fourier Frequency derivative (# of bins     */
  /*       the signal smears over during the observation).             */
  /*    'w' is the Fourier Frequency 2nd derivative (change in the     */
  /*       Fourier f-dot during the observation).                      */
  /*    'accuracy' is either LOWACC or HIGHACC.                        */
  /*  Notes:                                                           */
  /*    The result must be multiplied by 2*'numbetween' to get the     */
  /*    length of the array required to hold such a kernel.            */
    if (fabs(w) < 1.0e-7)
        return presto_z_resp_halfwidth(z, accuracy);
    double r0 = 0.5 * (w / 6.0 - z); // Starting deviation from r_avg
    double r1 = 0.5 * (w / 6.0 + z); // Ending deviation from r_avg
    // We need to know the maximum deviation from r_avg
    double maxdev = fabs(r0) > fabs(r1) ? fabs(r0) : fabs(r1);
    // If the extrema of the parabola is within 0 < u < 1, then
    // it will be a new freq minimum or maximum
    double u_ext = 0.5 - z / w;
    if (u_ext > 0.0 && u_ext < 1.0) {
        double z0 = z - w / 2.0; // Starting z
        // Value of r at the extremum
        double r_ext =  0.5 * w * u_ext * u_ext + z0 * u_ext + r0;
        maxdev = fabs(r_ext) > maxdev ? fabs(r_ext) : maxdev;
    }
    if (accuracy == HIGHACC) {
        return (int) (1.1 * maxdev) + NUMFINTBINS * 3;
    } else {
		return (int) (1.1 * maxdev) + NUMFINTBINS;
    }
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

    response = (cufftComplex*)malloc(numkern*sizeof(cufftComplex));
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


  cufftComplex* presto_gen_z_response(double roffset, int numbetween, double z, int numkern)
  /*  Generate the response function for Fourier f-dot interpolation.  */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 2 = interbins) */
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

    /* Check that the arguments are OK */

    if (roffset < 0.0 || roffset >= 1.0) {
        printf("\n  roffset = %f (out of bounds) in gen_z_response().\n\n", roffset);
		exit(-1);
	}
	if (numbetween < 1 || numbetween >= 20000) {
		printf("\n  numbetween = %d (out of bounds) in gen_z_response().\n\n", numbetween);
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

    /* If z~=0 use the normal Fourier interpolation kernel */
   
    absz = fabs(z);
    if (absz < 1E-4) {
     response = presto_gen_r_response(roffset, numbetween, numkern);
      return response;
    }
    response = (cufftComplex*)malloc(numkern*sizeof(cufftComplex));

    /* Begin the calculations */

    startr = roffset - (0.5 * z);
    startroffset = (startr < 0) ? 1.0 + modf(startr, &tmprl) : modf(startr, &tmprl);
    signz = (z < 0.0) ? -1 : 1;
    zd = signz * M_SQRT2 / sqrt(absz);
    cons = zd / 2.0;
    pibyz = M_PI / z;
    startr += numkern / (double)(2.0*numbetween);
    delta = -1.0 / numbetween;

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

    return response;
  }

// Tests conducted checking the fractional deviation of the amplitudes
// of the w-response calculation using different num_pts_wdat,
// compared to 262144.  roffset=[0,1], z=[-200,200], w=[-1000,1000]
//
// NUM_PTS_WDAT  MinFracDev   MedFracDev  MaxFracDev
//   131072      1.5983e-05   6.4267e-05   0.002060
//    65536      5.1875e-05   0.00021747   0.005147
//    32768      0.00012699   0.00051079   0.012568
//    16384      0.00027375   0.00112215   0.026279
//     8192      0.00054102   0.00221496   0.053507
//     4096      0.00104040   0.00410371   0.101785
//     2048      0.00244757   0.00875644   0.224530
//     1024      0.00427585   0.01669957   0.497524

cufftComplex *presto_gen_w_response(double roffset, int numbetween, double z, double w, int numkern){
  /*  Generate the response for Fourier f, f-dot, f-dotdot interp.     */
  /*  Arguments:                                                       */
  /*    'roffset' is the offset in Fourier bins for the full response  */
  /*       (i.e. At this point, the response would equal 1.0)          */
  /*    'numbetween' is the number of points to interpolate between    */
  /*       each standard FFT bin.  (i.e. 'numbetween' = 2 = interbins) */
  /*    'z' is the average Fourier Frequency derivative (# of bins     */
  /*       the signal smears over during the observation).             */
  /*    'w' is the Fourier Frequency 2nd derivative (change in the     */
  /*       Fourier f-dot during the observation).                      */
  /*    'numkern' is the number of complex points that the kernel will */
  /*       contain.                                                    */
  /*  This version uses zero-padding to get the "numbetween"           */
    int ii, fbar, num_pts_wdat;
    float *data;
	cufftComplex *data_fft;
    double amp, f, fd, fdd, dt, t, phase, dfbar;
    cufftComplex *response;

    /* Check that the arguments are OK */
    if (roffset < 0.0 || roffset >= 1.0) {
        printf("\n  roffset = %f (out of bounds) in gen_w_response().\n\n", roffset);
        exit(-1);
    }
    if (numbetween < 1 || numbetween >= 20000) {
        printf("\n  numbetween = %d (out of bounds) in gen_w_response().\n\n",
               numbetween);
        exit(-1);
    }
    if (numkern < numbetween) {
        printf("\n  numkern = %d (out of bounds) in gen_w_response().\n\n", numkern);
        exit(-1);
    }
    if ((numkern % (2 * numbetween)) != 0) {
        printf("\n  numkern %% (2 * numbetween) != 0 in gen_w_response().\n\n");
        exit(-1);
    }

    /* If w~=0 use the normal F-dot Fourier interpolation kernel */
    if (fabs(w) < 1E-4) {
        response = presto_gen_z_response(roffset, numbetween, z, numkern);
        return response;
    }

    /* Choose num_pts_wdat so that there is plenty of Freq range */
    /* outside of the RZW response. */
    num_pts_wdat = next2_to_n(6 * presto_w_resp_halfwidth(z, w, LOWACC) + 200 + numkern / numbetween);

    /* Otherwise initialize some data */
    dt = 1.0 / (double) num_pts_wdat;
    amp = 2.0 * dt;
    fbar = num_pts_wdat / 4;  // num_pts_wdat / 4 is average freq
    dfbar = (double) fbar + roffset;
    // r_o = rbar - zbar/2 + w/12  where _o is initial and bar is average
    // z_o = zbar - w/2
    f = dfbar - 0.5 * z + w / 12.0;     //  This shifts the initial f appropriately
    fd = (z - 0.5 * w) / 2.0;   // z - w/2 is the initial z value
    fdd = w / 6.0;

	/* Generate the data set.  Use zero-padding to do the interpolation. */
	size_t data_size = num_pts_wdat*numbetween*sizeof(float);
	size_t data_fft_size = num_pts_wdat*numbetween*sizeof(cufftComplex);
	data = (float*) malloc(data_size);
	if(data==NULL) printf("Error allocating data!\n");
	data_fft = (cufftComplex*) malloc(data_fft_size);
	if(data_fft==NULL) printf("Error allocating data_fft!\n");
	
    for (ii = 0; ii < num_pts_wdat * numbetween; ii++){
        data[ii] = 0.0;
	}
    for (ii = 0; ii < num_pts_wdat; ii++) {
        t = ii * dt;
        phase = TWOPI * (t * (t * (t * fdd + fd) + f));
        data[ii] = amp * cos(phase);
    }
    
	/* FFT the data */
	//---------- CUFFT ------------->
	cufftHandle plan;
	cufftResult cuFFT_error;
	cudaError_t cudaError;
	cuFFT_error = cufftPlan1d(&plan, num_pts_wdat*numbetween, CUFFT_R2C, 1);
	if (CUFFT_SUCCESS == cuFFT_error) {
		float *d_input;
		cudaMalloc((void **) &d_input, data_size);
		cufftComplex *d_output;
		cudaMalloc((void **) &d_output, data_fft_size);
		
		cudaError = cudaMemcpy(d_input, data, data_size, cudaMemcpyHostToDevice);
	    if(cudaError != cudaSuccess) printf("Could not cudaMemcpy in presto_func.cpp");
	    
		cufftExecR2C(plan, d_input, d_output);
		
		cudaError = cudaMemcpy(data_fft, d_output, data_fft_size, cudaMemcpyDeviceToHost);
	    if(cudaError != cudaSuccess) printf("Could not cudaMemcpy in presto_func.cpp");
		
		cudaFree(d_input);
		cudaFree(d_output);
	}
	else printf("CUFFT error: Plan creation failed");
	cufftDestroy(plan);
	//------------------------------<
	
	//--------- fftw --------------->
    //fftwf_plan rffplan = fftwf_plan_dft_r2c_1d(num_pts_wdat*numbetween, data, (fftwf_complex*) data_fft, FFTW_ESTIMATE);
	//if(rffplan==NULL) printf("FFT error!\n");
	//else {
	//	fftwf_execute(rffplan);
	//	fftwf_destroy_plan(rffplan);
	//}
	//------------------------------<
	
    /* Generate the final response */
    response = (cufftComplex*) malloc(numkern * sizeof(cufftComplex));
	if(response==NULL) printf("Error allocating response!\n");
    /* Chop off the contaminated ends and/or the extra data */
    memcpy(response, data_fft + (fbar * numbetween - numkern / 2), sizeof(cufftComplex)*numkern);
	
	
    /* cleanup */
    free(data);
	free(data_fft);
    return response;
}

void presto_place_complex_kernel(cufftComplex * kernel, int numkernel, cufftComplex * result, int numresult)
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
      result[ii] = zeros;
    }

    memcpy(result, kernel + halfwidth, sizeof(cufftComplex) * halfwidth);
    memcpy(result + numresult - halfwidth, kernel, sizeof(cufftComplex) * halfwidth);
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

} //namespace astroaccelerate
