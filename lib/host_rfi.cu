#include <stdio.h>
#include <stdlib.h>
#include "headers/params.h"

void rfi(int nsamp, int nchans, unsigned short **input_buffer)
{
/*
	int chan_counter = 0;
	int good = -1;
	int flag = -1;

	 for(int t=0; t<nsamp; t++) {
	 for(int c = 0; c < nchans-30; c++) {
	 float running_mean=0.0f;
	 for(int w=0; w<30; w++) {
	 running_mean+=(*input_buffer)[c + w + nchans*(t)];
	 }
	 running_mean=running_mean/30.0f;
	 float sd_sum=0.0f;
	 for(int w=0; w<30; w++) {
	 sd_sum+=(((*input_buffer)[c + w+ nchans*(t)]-running_mean)*((*input_buffer)[c + w+ nchans*(t)]-running_mean));
	 }
	 sd_sum=sqrt(sd_sum/30.0f);
	 float test_mean=0.0f;
	 for(int w=0; w<5; w++) {
	 test_mean+=(*input_buffer)[c + w+ nchans*(t)];
	 }
	 test_mean=test_mean/5.0f;
	 if((test_mean-running_mean) > 2.5*sd_sum) {
	 if(good==-1) {
	 for(int j=0; j<nchans; j++) {
	 (*input_buffer)[j + nchans*t]=127;
	 }
	 } else {
	 for(int j=0; j<nchans; j++) {
	 (*input_buffer)[j + nchans*t]=(*input_buffer)[good*nchans +(int)((rand()/((float)RAND_MAX))*nchans)];
	 }
	 }
	 printf("\nclipping:\t%d, %d, %d\t %f %f %f", c, good, t, test_mean, running_mean, sd_sum);
	 chan_counter++;
	 flag=0;
	 break;
	 }
	 }
	 if(flag == -1) good = t;
	 flag = -1;
	 }
	 printf("\nClipped %lf percent of spectra", ((double)chan_counter/(double)nchans)*100.0);

	 chan_counter=0;
	 good=-1;
	 flag=-1;
	 for(int c = 0; c < nchans; c++) {
	 for(int t=0; t<nsamp-300; t++) {
	 float running_mean=0.0f;
	 for(int w=0; w<300; w++) {
	 running_mean+=(*input_buffer)[c + nchans*(t + w)];
	 }
	 running_mean=running_mean/300.0f;
	 float sd_sum=0.0f;
	 for(int w=0; w<300; w++) {
	 sd_sum+=(((*input_buffer)[c + nchans*(t + w)]-running_mean)*((*input_buffer)[c + nchans*(t + w)]-running_mean));
	 }
	 sd_sum=sqrt(sd_sum/300.0f);
	 float test_mean=0.0f;
	 for(int w=0; w<5; w++) {
	 test_mean+=(*input_buffer)[c + nchans*(t + w)];
	 }
	 test_mean=test_mean/5.0f;
	 if((test_mean-running_mean) > 4.25*sd_sum) {
	 if(good==-1) {
	 for(int j=0; j<nsamp; j++) {
	 (*input_buffer)[c + nchans*j]=127.3858f;
	 }
	 } else {
	 for(int j=0; j<nsamp; j++) {
	 (*input_buffer)[c + nchans*j]=(*input_buffer)[good + (int)((rand()/((float)RAND_MAX))*(nsamp-1))*nchans];
	 //printf("\n%d", (int)((rand()/(float)RAND_MAX)*(nsamp-1)));
	 //(*input_buffer)[c*nsamp + j]=127.3858f;
	 }
	 }
	 printf("\nclipping:\t%d, %d, %d\t %f %f %f", c, good, t, test_mean, running_mean, sd_sum);
	 chan_counter++;
	 flag=0;
	 break;
	 }
	 }
	 if(flag == -1) good = c;
	 flag = -1;
	 }
	 printf("\nClipped %lf percent of channels", ((double)chan_counter/(double)nchans)*100.0);
	 */
	// Zero dm 
//#pragma omp parallel for
	for (int t = 0; t < nsamp; t++)
	{
		double mean = 0.0;
		for (int c = 0; c < nchans; c++)
		{
			mean += (double) ( *input_buffer )[c + (unsigned long) nchans * t];
		}
		mean = mean / (double) nchans;
		for (int c = 0; c < nchans; c++)
		{
			( *input_buffer )[c + (unsigned long) nchans * t] = (unsigned short) ( (unsigned char) ( (double) ( *input_buffer )[c + (unsigned long) nchans * t] - mean ) );
		}
	}

	/*

	 int chan_counter=0;
	 int good=-1;
	 int flag=-1;
	 for(int c = 0; c < nchans; c++) {
	 total  = 0.0;
	 counter = 0;
	 double sd_total=0.0;
	 double sd=0.0;
	 for(int t=0; t<nsamp; t++) {
	 total += (double)(*input_buffer)[c*nsamp + t];
	 counter++;
	 }
	 mean = ((double)total)/((double)counter);  // Mean for data sample
	 for(int t=0; t<nsamp; t++) {
	 total = ((double)(*input_buffer)[c*nsamp + t]-mean);
	 sd_total += (total*total);
	 //		(*input_buffer)[c*nsamp + t] = (float)(total);
	 }
	 sd=(sqrt(sd_total / (double)nsamp));
	 //printf("\nc:\t%d, m:\t%lf, sd:\t%lf", c, mean, sd);
	 float clip_factor=(254.0f-mean)/sd;

	 */
	/*
	 for(int t=0; t<nsamp-300; t++) {
	 //if(((*input_buffer)[c*nsamp + t]-mean) > clip_factor*sd &&(*input_buffer)[c*nsamp + t] != 0.0) {
	 float running_mean=0.0f;
	 for(int w=0; w<300; w++) {
	 running_mean+=(*input_buffer)[c*nsamp + t + w];
	 }
	 running_mean=running_mean/300.0f;
	 float sd_sum=0.0f;
	 for(int w=0; w<300; w++) {
	 sd_sum+=(((*input_buffer)[c*nsamp + t + w]-running_mean)*((*input_buffer)[c*nsamp + t + w]-running_mean));
	 }
	 sd_sum=sqrt(sd_sum/300.0f);
	 float test_mean=0.0f;
	 for(int w=0; w<5; w++) {
	 test_mean+=(*input_buffer)[c*nsamp + t + w];
	 }
	 test_mean=test_mean/5.0f;
	 if((test_mean-running_mean) > 4.5*sd_sum) {
	 if(good==-1) {
	 for(int j=0; j<nsamp; j++) {
	 (*input_buffer)[c*nsamp + j]=127.3858f;
	 }
	 } else {
	 for(int j=0; j<nsamp; j++) {
	 (*input_buffer)[c*nsamp + j]=(*input_buffer)[good*nsamp + (int)((rand()/((float)RAND_MAX))*(nsamp-1))];
	 //printf("\n%d", (int)((rand()/(float)RAND_MAX)*(nsamp-1)));
	 //(*input_buffer)[c*nsamp + j]=127.3858f;
	 }
	 }
	 printf("\nclipping:\t%d, %d, %d\t %f %f %f", c, good, t, test_mean, running_mean, sd_sum);
	 chan_counter++;
	 flag=0;
	 break;
	 }
	 }
	 if(flag == -1) good = c;
	 flag = -1;
	 }
	 printf("\nClipped %lf percent of channels", ((double)chan_counter/(double)nchans)*100.0);
	 */
	/*
	 for(int t=0; t<nsamp; t++) {
	 total  = 0.0;
	 counter = 0;
	 for(int c = 0; c < nchans; c++) {

	 total += (*input_buffer)[t*nchans + c];
	 counter++;
	 }
	 mean = ((double)total)/((double)counter);  // Mean for data sample
	 for(int c = 0; c < nchans; c++) {
	 (*input_buffer)[t*nchans + c] = (float)(((double)(*input_buffer)[t*nchans + c]-mean));
	 }
	 }
	 */
	/*
	 unsigned long int j;
	 unsigned long int vals;
	 unsigned long int counter;
	 unsigned long int mean_counter;

	 double mean, mean_sub,  stddev, mean_per_channel;

	 double	total;

	 // Calculate the total number of values
	 vals = (unsigned long int)(nsamp*nchans);

	 // Calculate the mean
	 total  = 0.0;
	 #pragma omp parallel for default(shared) private(j) reduction(+:total)
	 for(j = 0; j < vals; j++) {
	 total += (double)(*input_buffer)[j];
	 }
	 mean = (total/(double)vals);  // Mean for data sample

	 //	// Calculate standard deviation
	 //	total = 0.0;
	 //	#pragma omp parallel for default(shared) private(j) reduction(+:total)
	 //	for(j = 0; j < vals; j++) {
	 //		total += (double)(((double)(*input_buffer)[j] - mean)*((double)(*input_buffer)[j] - mean));
	 //	}
	 //	stddev = sqrt(total / (double)vals); // Stddev for data sample
	 //
	 //	printf("\nMean: %lf, Stddev: %lf", mean, stddev), fflush(stdout);

	 // Calculate the mean per channel
	 //#pragma omp parallel for default(shared) private(j) reduction(+:total)
	 double sd_total=0.0;
	 for(int c = 0; c < nchans; c++) {
	 total  = 0.0;
	 counter = 0;
	 for(int t=0; t<nsamp; t++) {
	 total += (double)(*input_buffer)[c*nsamp+t];
	 counter++;
	 }
	 mean_per_channel = (total/(double)counter);  // Mean for data sample

	 sd_total += (double)(((double)(mean_per_channel - mean)*((double)(mean_per_channel - mean))));
	 }
	 stddev = sqrt(sd_total / (double)nchans); // Stddev for data sample

	 mean_counter=0;
	 for(int c = 0; c < nchans; c++) {
	 total  = 0.0;
	 counter = 0;
	 for(int t=0; t<nsamp; t++) {
	 total += (double)(*input_buffer)[c*nsamp+t];
	 counter++;
	 }
	 mean_per_channel = (total/(double)counter);  // Mean for data sample

	 if(abs(mean_per_channel - mean) < 2*stddev) {
	 mean_sub+= mean_per_channel;
	 mean_counter++;
	 }
	 }
	 mean = (mean_sub/(double)mean_counter);  // Mean for data sample

	 for(int c = 0; c < nchans; c++) {
	 total  = 0.0;
	 counter = 0;
	 for(int t=0; t<nsamp; t++) {
	 total += (double)(*input_buffer)[c*nsamp+t];
	 counter++;
	 }
	 mean_per_channel = (total/(double)counter);  // Mean for data sample

	 if(abs(mean_per_channel - mean) > 2*stddev) {
	 printf("\n Striking out channnel:\t%d", c);
	 for(int t=0; t<nsamp; t++) {
	 (*input_buffer)[c*nsamp+t] = (unsigned short)mean;
	 }
	 }
	 }
	 */
	/*
	 // Calculate the total number of values
	 vals = (unsigned long int)(nsamp*nchans);

	 // Calculate the mean
	 total  = 0.0;
	 #pragma omp parallel for default(shared) private(j) reduction(+:total)
	 for(j = 0; j < vals; j++) {
	 total += (double)(*input_buffer)[j];
	 }
	 mean = (total/(double)vals);  // Mean for data sample

	 // Calculate standard deviation
	 total = 0.0;
	 #pragma omp parallel for default(shared) private(j) reduction(+:total)
	 for(j = 0; j < vals; j++) {
	 total += (double)(((double)(*input_buffer)[j] - mean)*((double)(*input_buffer)[j] - mean));
	 }
	 stddev_orig = sqrt(total / (double)vals); // Stddev for data sample

	 printf("\nMean: %lf, Stddev: %lf", mean, stddev_orig), fflush(stdout);


	 unsigned short sd_h=(unsigned short)10*stddev_orig;
	 unsigned short sd_l=(unsigned short)4*stddev_orig;

	 #pragma omp parallel for
	 for(j = 0; j < vals; j++) {
	 if((*input_buffer)[j] > sd_h) (*input_buffer)[j]=mean;
	 //if((*input_buffer)[j] < sd_l) (*input_buffer)[j]=mean;
	 }
	 */
}
