# **Single Pulse Search**

Single pulse search consists from several stages, currently these are: calculation of standard deviation, SPDT and peak-finding. The description of how is mean and standard deviation calculated by AA is [here](calculation_of_mean_and_standard_deviation). 

We have implemented single pulse search algorithm which is rely on combination of boxcar filters of different widths and decimation of time-domain data in time. The algorithm we have used is configurable and could be adjusted to for speed or sensitivity. Sensitivity for currently used configuration is shown below. The figure shows two limits for the algorithm: R-SNRmax is best possible SNR recovered; R-SNRmin worst possible SNR recovered. 

![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/BOXDIT_R-SNR_32_16_8_integrated.jpg)

### Results
Results for signals generated by sigproc fake with different periods, DMs and widths are shown below.

![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/SPDT_P70_DM70_W43_SNR8.jpg)
Periodic with P=70ms, DM=70, initial SNR=8 and width w=43,75 samples (4% duty cycle) as detected by AA. Top right: detected pulses DM/time[s]; Top left: detected pulse width DM/width[samples]; Bottom right: time[s]/R-SNR; Bottom left: individual pulses DM/time/SNR;


![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/SPDT_P70_DM2100_W43_SNR8.jpg)
Periodic with P=70ms, DM=2100, initial SNR=8 and width w=43,75 samples (4% duty cycle) as detected by AA. Top right: detected pulses DM/time[s]; Top left: detected pulse width DM/width[samples]; Bottom right: time[s]/R-SNR; Bottom left: individual pulses DM/time/SNR;


![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/figures/SPDT_Wvar.jpg)
Reported pulse width for different signals at DM=90 with initial SNR=8. Top left: w=100samples; Top right: w=200 samples; Bottom right: w=400 samples; Bottom left: w=800 samples.

![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/SPDT_P100_DM50_W62_SNRvar.jpg)
Detection of periodic signals with different initial SNR (DM=50; P=100ms). Top left: SNR=8; Top right: SNR=4; Bottom right: SNR=2; Bottom left: SNR=1.