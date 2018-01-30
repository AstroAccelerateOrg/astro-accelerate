#include <stdio.h>
#include <stdlib.h>
#include "headers/params.h"

void rfi(int nsamp, int nchans, unsigned short **input_buffer)
{

	int	file_reducer 	    = 10;

	// Need to be adaptive and automatic.
	int	outlier_itterations = 10;
	float 	sigma_cut 	    = 2.0f;
	
	float 	*stage = (float*)malloc((size_t)nsamp*(size_t)nchans*sizeof(float));

	for(int c = 0; c < nchans; c++) {
		for(int t = 0; t < (nsamp); t++) {
			stage[c * (size_t)nsamp + t] = (float) (*input_buffer)[c  + (size_t)nchans * t];
		}
	}


	// ~~~ Output the raw data ~~~ //
	//FILE *fp_raw = fopen ("raw_chans.txt", "w+");
	//for(int c = 0; c < nchans; c++) {
	//	for(int t = 0; t < (nsamp)/file_reducer; t++) {
	//		fprintf(fp_raw, "%f ", stage[c * (size_t)nsamp + t]);
	//	}
	//	fprintf(fp_raw, "\n");
	//}
   	//fclose(fp_raw);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //


	
	// ~~~ RFI Correct ~~~ //	

	double orig_mean = 0.0;
	double orig_var=0.0;

        for(int c = 0; c < nchans; c++) {
       	        for(int t = 0; t < (nsamp); t++) orig_mean+=stage[c * (size_t)nsamp + t];
       	}
	orig_mean/=(nsamp*nchans);

        for(int c = 0; c < nchans; c++) {
		for(int t = 0; t < (nsamp); t++) {
			orig_var+=(stage[c * (size_t)nsamp + t]-orig_mean)*(stage[c * (size_t)nsamp + t]-orig_mean);
		}
	}
	orig_var/=(nsamp*nchans);
	orig_var=sqrt(orig_var);
	
	int *mask = (int*)malloc(nchans*sizeof(int));
	for(int c = 0; c < nchans; c++) mask[c]=1;

	int *spectra_mask = (int*)malloc(nsamp*sizeof(int));
	for(int t = 0; t < nsamp; t++) spectra_mask[t]=1;
 
	int *chan_mask_weights = (int*)malloc(nchans*sizeof(int));
	for(int c = 0; c < nchans; c++) chan_mask_weights[c]=0;

	int *spectra_mask_weights = (int*)malloc(nsamp*sizeof(int));
	for(int t = 0; t < nsamp; t++) spectra_mask_weights[t]=0;

	double *chan_mean = (double*)malloc(nchans*sizeof(double));
       	for(int c = 0; c < nchans; c++) chan_mean[c] = 0.0;

	double *chan_var = (double*)malloc(nsamp*sizeof(double));
       	for(int c = 0; c < nchans; c++) chan_var[c] = 0.0;

	double *spectra_mean = (double*)malloc(nsamp*sizeof(double));
       	for(int t = 0; t < nsamp; t++) spectra_mean[t] = 0.0;

	double *spectra_var = (double*)malloc(nsamp*sizeof(double));
       	for(int t = 0; t < nsamp; t++) spectra_var[t] = 0.0;

        for(int c = 0; c < nchans; c++) {
  
		int counter = 0;

		for(int t = 0; t < nsamp; t++) spectra_mask[t]=1;

		for(int r = 0; r < outlier_itterations; r++) {
		
			counter = 0;
			chan_mean[c] = 0.0;
        	        for(int t = 0; t < (nsamp); t++) {
				if(spectra_mask[t] == 1) {
	                	        chan_mean[c]+=stage[c * (size_t)nsamp + t];
					counter++;
				}
        	        }
			if(counter == 0) {
				//printf("\nCounter zero, Channel %d", c);
				break;
			}
			chan_mean[c]/=(counter);



			counter=0;
			chan_var[c]=0.0;
			for(int t = 0; t < (nsamp); t++) {
				if(spectra_mask[t] == 1) {
					chan_var[c]+=(stage[c * (size_t)nsamp + t]-chan_mean[c])*(stage[c * (size_t)nsamp + t]-chan_mean[c]);
					counter++;
				}
			}
			if(counter == 0) {
				//printf("\nCounter zero, Channel %d", c);
				break;
			}
			chan_var[c]/=(counter);
			chan_var[c]=sqrt(chan_var[c]);

			if(abs(chan_var[c])*1000000.0 < 0.1) {
				//printf("\nVarience zero, Channel %d", c);
				chan_var[c] = 1.0;
				break;
			}

			for(int t = 0; t < (nsamp); t++) {
				if(((stage[c * (size_t)nsamp + t]-chan_mean[c])/chan_var[c]) > sigma_cut || ((stage[c * (size_t)nsamp + t]-chan_mean[c])/chan_var[c]) < -sigma_cut) {
					spectra_mask[t]=0;
				} else {
					spectra_mask[t]=1;
					spectra_mask_weights[t]++;
				}
			}
		}

		if(counter != 0) {
			for(int t = 0; t < (nsamp); t++) {
				stage[c * (size_t)nsamp + t]=(stage[c * (size_t)nsamp + t]-(float)chan_mean[c])/(float)chan_var[c];
			}
		}
	}

	for(int t = 0; t < (nsamp); t++) {

		int counter = 0;

		for(int c = 0; c < nchans; c++) mask[c]=1;

		for(int r = 0; r < outlier_itterations; r++) {
		
			counter = 0;
			spectra_mean[t] = 0.0;
			for(int c = 0; c < nchans; c++) {
				if(mask[c] == 1) {
	                	        spectra_mean[t]+=stage[c * (size_t)nsamp + t];
					counter++;
				}
        	        }
			if(counter == 0) {
				//printf("\nCounter zero, Spectra %d", t);
				break;
			}
			spectra_mean[t]/=(counter);

			counter=0;
			spectra_var[t]=0.0;
			for(int c = 0; c < nchans; c++) {
				if(mask[c] == 1) {
					spectra_var[t]+=(stage[c * (size_t)nsamp + t]-spectra_mean[t])*(stage[c * (size_t)nsamp + t]-spectra_mean[t]);
					counter++;
				}
			}
			if(counter == 0) {
				//printf("\nCounter zero, Spectra %d", t);
				break;
			}
			spectra_var[t]/=((double)(counter));
			spectra_var[t]=sqrt(spectra_var[t]);


			if(abs(spectra_var[t])*1000000.0 < 0.1) {
				//printf("\nVarience zero, Spectra %d", t);
				spectra_var[t] = 1.0;
				break;
			}

			for(int c = 0; c < nchans; c++) {
				if(((stage[c * (size_t)nsamp + t]-spectra_mean[t])/spectra_var[t]) > sigma_cut || ((stage[c * (size_t)nsamp + t]-spectra_mean[t])/spectra_var[t]) < -sigma_cut) {
					mask[c]=0;
				} else {
					mask[c]=1;
					chan_mask_weights[c]++;
				}
			}
		}

		if(counter != 0) {
			for(int c = 0; c < nchans; c++) {
				stage[c * (size_t)nsamp + t]=(stage[c * (size_t)nsamp + t]-(float)spectra_mean[t])/(float)spectra_var[t];
			}
		}
	}


	int best_weight = 0;
	int best_spectra = 0;	
	for(int t = 0; t < (nsamp); t++) {
		if(best_weight < spectra_mask_weights[t]) best_spectra = t;
	}

	best_weight = 0;
	int best_chan = 0;	
	for(int c = 0; c < nchans; c++) {
		if(best_weight < chan_mask_weights[c]) best_chan = c;
	}

        for(int c = 0; c < nchans; c++) {
		//if(mask[c] == 0 || chan_var[c] == 1.0) {
		if(mask[c] == 0) {
			int perm_one = (int)(((float)rand()/(float)RAND_MAX)*nsamp)+1;
			for(int t = 0; t < (nsamp); t++) {
				stage[(c) * (size_t)nsamp + t] = stage[(best_chan) * (size_t)nsamp + (t+perm_one)%nsamp];
			}
		}
	}

	for(int t = 0; t < (nsamp); t++) {
		//if(spectra_mask[t] == 0 || spectra_var[t] == 1.0) {
		if(spectra_mask[t] == 0) {
			int perm_one = (int)(((float)rand()/(float)RAND_MAX)*nchans)+1;
			for(int c = 0; c < nchans; c++) {
				stage[(c) * (size_t)nsamp + t] = stage[(c+perm_one)%nchans * (size_t)nsamp + (best_spectra)];
			}
		}
	}


	// Global mask
	
	double mean = 0.0;
	double var = 0.0;
	for(int t = 0; t < (nsamp); t++) {
		mean+=spectra_mean[t];
	}
	mean/=nsamp;

	for(int t = 0; t < (nsamp); t++) {
		var+=(spectra_mean[t] - mean)*(spectra_mean[t] - mean);
	}
	var/=(nsamp);
	var=sqrt(var);
	
	int counter = 0;
	for(int t = 0; t < (nsamp); t++) {
		if(((spectra_mean[t] - mean)/var) > sigma_cut || ((spectra_mean[t] - mean)/var) < -sigma_cut) {
			spectra_mask[t]=0;
			//printf("\nmasking spectra:\t%d", t);
			counter++;
		} else {
			spectra_mask[t]=1;
		}
	}
	printf("\nspectra masked:\t%f", ((float)counter/(float)nsamp)*100);

	for(int t = 0; t < (nsamp); t++) {
		if(spectra_mask[t] == 0) {
			int perm_one = (int)(((float)rand()/(float)RAND_MAX)*nchans)+1;
			for(int c = 0; c < nchans; c++) {
				stage[(c) * (size_t)nsamp + t] = stage[(c+perm_one)%nchans * (size_t)nsamp + (best_spectra)];
			}
		}
	}


	mean = 0.0;
	var = 0.0;
	for(int c = 0; c < nchans; c++) {
		mean+=chan_mean[c];
	}
	mean/=nchans;

	for(int c = 0; c < nchans; c++) {
		var+=(chan_mean[c] - mean)*(chan_mean[c] - mean);
	}
	var/=(nchans);
	var=sqrt(var);

	counter = 0;
	for(int c = 0; c < nchans; c++) {
		if(((chan_mean[c] - mean)/var) > sigma_cut || ((chan_mean[c] - mean)/var) < -sigma_cut) {
			mask[c]=0;
			//printf("\nmasking chan:\t%d", c);
			counter++;
		} else {
			mask[c]=1;
		}
	}
	printf("\nchannels masked:\t%f", ((float)counter/(float)nchans)*100);

        for(int c = 0; c < nchans; c++) {
		if(mask[c] == 0) {
			int perm_one = (int)(((float)rand()/(float)RAND_MAX)*nsamp)+1;
			for(int t = 0; t < (nsamp); t++) {
				stage[(c) * (size_t)nsamp + t] = stage[(best_chan) * (size_t)nsamp + (t+perm_one)%nsamp];
			}
		}
	}


	for(int c = 0; c < nchans; c++) {
		for(int t = 0; t < (nsamp); t++) {
			(*input_buffer)[c  + (size_t)nchans * t] = (unsigned char) ((stage[c * (size_t)nsamp + t]*orig_var)+orig_mean);
		}
	}



	//FILE *fp_mask = fopen ("masked_chans.txt", "w+");
	//for(int c = 0; c < nchans-1; c++) {
	//	for(int t = 0; t < (nsamp)/file_reducer; t++) {
	//		fprintf(fp_mask, "%f ", ((stage[c * (size_t)nsamp + t]*25.6)+128));
	//	}
	//	fprintf(fp_mask, "\n");
	//}
   	//fclose(fp_mask);



	free(mask);
	free(spectra_mask);
	free(chan_mask_weights);
	free(spectra_mask_weights);
	free(chan_mean);
	free(chan_var);
	free(spectra_mean);
	free(spectra_var);
	free(stage);
}

	
	/*
	int chan_counter = 0;
	int good = -1;
	int flag = -1;

	 for(int t=0; t<nsamp; t++) {
	 for(int c = 0; c < nchans-30; c++) {
	 float running_mean=0.0f;
	 for(int w=0; w<30; w++) {
	 running_mean+=stage[c + w + nchans*(t)];
	 }
	 running_mean=running_mean/30.0f;
	 float sd_sum=0.0f;
	 for(int w=0; w<30; w++) {
	 sd_sum+=((stage[c + w+ nchans*(t)]-running_mean)*(stage[c + w+ nchans*(t)]-running_mean));
	 }
	 sd_sum=sqrt(sd_sum/30.0f);
	 float test_mean=0.0f;
	 for(int w=0; w<5; w++) {
	 test_mean+=stage[c + w+ nchans*(t)];
	 }
	 test_mean=test_mean/5.0f;
	 if((test_mean-running_mean) > 2.5*sd_sum) {
	 if(good==-1) {
	 for(int j=0; j<nchans; j++) {
	 stage[j + nchans*t]=127;
	 }
	 } else {
	 for(int j=0; j<nchans; j++) {
	 stage[j + nchans*t]=stage[good*nchans +(int)((rand()/((float)RAND_MAX))*nchans)];
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
	 running_mean+=stage[c + nchans*(t + w)];
	 }
	 running_mean=running_mean/300.0f;
	 float sd_sum=0.0f;
	 for(int w=0; w<300; w++) {
	 sd_sum+=((stage[c + nchans*(t + w)]-running_mean)*(stage[c + nchans*(t + w)]-running_mean));
	 }
	 sd_sum=sqrt(sd_sum/300.0f);
	 float test_mean=0.0f;
	 for(int w=0; w<5; w++) {
	 test_mean+=stage[c + nchans*(t + w)];
	 }
	 test_mean=test_mean/5.0f;
	 if((test_mean-running_mean) > 4.25*sd_sum) {
	 if(good==-1) {
	 for(int j=0; j<nsamp; j++) {
	 stage[c + nchans*j]=127.3858f;
	 }
	 } else {
	 for(int j=0; j<nsamp; j++) {
	 stage[c + nchans*j]=stage[good + (int)((rand()/((float)RAND_MAX))*(nsamp-1))*nchans];
	 //printf("\n%d", (int)((rand()/(float)RAND_MAX)*(nsamp-1)));
	 //stage[c*nsamp + j]=127.3858f;
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
	//for (int t = 0; t < nsamp; t++)
	//{
	//	double mean = 0.0;
	//	for (int c = 0; c < nchans; c++)
	//	{
	//		mean += (double) ( *input_buffer )[c + (unsigned long) nchans * t];
	//	}
	//	mean = mean / (double) nchans;
	//	for (int c = 0; c < nchans; c++)
	//	{
	//		( *input_buffer )[c + (unsigned long) nchans * t] = (unsigned short) ( (unsigned char) ( (double) ( *input_buffer )[c + (unsigned long) nchans * t] - mean ) );
	//	}
	//}

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
	 total += (double)stage[c*nsamp + t];
	 counter++;
	 }
	 mean = ((double)total)/((double)counter);  // Mean for data sample
	 for(int t=0; t<nsamp; t++) {
	 total = ((double)stage[c*nsamp + t]-mean);
	 sd_total += (total*total);
	 //		stage[c*nsamp + t] = (float)(total);
	 }
	 sd=(sqrt(sd_total / (double)nsamp));
	 //printf("\nc:\t%d, m:\t%lf, sd:\t%lf", c, mean, sd);
	 float clip_factor=(254.0f-mean)/sd;

	 */
	/*
	 for(int t=0; t<nsamp-300; t++) {
	 //if((stage[c*nsamp + t]-mean) > clip_factor*sd &&stage[c*nsamp + t] != 0.0) {
	 float running_mean=0.0f;
	 for(int w=0; w<300; w++) {
	 running_mean+=stage[c*nsamp + t + w];
	 }
	 running_mean=running_mean/300.0f;
	 float sd_sum=0.0f;
	 for(int w=0; w<300; w++) {
	 sd_sum+=((stage[c*nsamp + t + w]-running_mean)*(stage[c*nsamp + t + w]-running_mean));
	 }
	 sd_sum=sqrt(sd_sum/300.0f);
	 float test_mean=0.0f;
	 for(int w=0; w<5; w++) {
	 test_mean+=stage[c*nsamp + t + w];
	 }
	 test_mean=test_mean/5.0f;
	 if((test_mean-running_mean) > 4.5*sd_sum) {
	 if(good==-1) {
	 for(int j=0; j<nsamp; j++) {
	 stage[c*nsamp + j]=127.3858f;
	 }
	 } else {
	 for(int j=0; j<nsamp; j++) {
	 stage[c*nsamp + j]=stage[good*nsamp + (int)((rand()/((float)RAND_MAX))*(nsamp-1))];
	 //printf("\n%d", (int)((rand()/(float)RAND_MAX)*(nsamp-1)));
	 //stage[c*nsamp + j]=127.3858f;
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

	 total += stage[t*nchans + c];
	 counter++;
	 }
	 mean = ((double)total)/((double)counter);  // Mean for data sample
	 for(int c = 0; c < nchans; c++) {
	 stage[t*nchans + c] = (float)(((double)stage[t*nchans + c]-mean));
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
	 total += (double)stage[j];
	 }
	 mean = (total/(double)vals);  // Mean for data sample

	 //	// Calculate standard deviation
	 //	total = 0.0;
	 //	#pragma omp parallel for default(shared) private(j) reduction(+:total)
	 //	for(j = 0; j < vals; j++) {
	 //		total += (double)(((double)stage[j] - mean)*((double)stage[j] - mean));
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
	 total += (double)stage[c*nsamp+t];
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
	 total += (double)stage[c*nsamp+t];
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
	 total += (double)stage[c*nsamp+t];
	 counter++;
	 }
	 mean_per_channel = (total/(double)counter);  // Mean for data sample

	 if(abs(mean_per_channel - mean) > 2*stddev) {
	 printf("\n Striking out channnel:\t%d", c);
	 for(int t=0; t<nsamp; t++) {
	 stage[c*nsamp+t] = (unsigned short)mean;
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
	 total += (double)stage[j];
	 }
	 mean = (total/(double)vals);  // Mean for data sample

	 // Calculate standard deviation
	 total = 0.0;
	 #pragma omp parallel for default(shared) private(j) reduction(+:total)
	 for(j = 0; j < vals; j++) {
	 total += (double)(((double)stage[j] - mean)*((double)stage[j] - mean));
	 }
	 stddev_orig = sqrt(total / (double)vals); // Stddev for data sample

	 printf("\nMean: %lf, Stddev: %lf", mean, stddev_orig), fflush(stdout);


	 unsigned short sd_h=(unsigned short)10*stddev_orig;
	 unsigned short sd_l=(unsigned short)4*stddev_orig;

	 #pragma omp parallel for
	 for(j = 0; j < vals; j++) {
	 if(stage[j] > sd_h) stage[j]=mean;
	 //if(stage[j] < sd_l) stage[j]=mean;
	 }
	 */
//}
