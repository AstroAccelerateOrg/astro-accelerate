# **Dedispersion kernel parameter optimisation**

It is possible to run the dedispersion step of AstroAccelerate using optimised values for each dedispersion measure (DM) range used, for the following parameters:
* SDIVINT - number of timesteps per block
* SDIVINDM - number of dedispersion measures per block
* UNROLLS
* SNUMREG

However, this was generally not found to be worthwhile due to adding complexity to the code for very small gains in performance.  

The optimised values are found by adjusting the dedispersion step to work only on one DM range at a time and record the speedup, then running the program multiple times for each range with a different set of parameters each time until the best speedup for each DM range is found. 

# Results
## Parameter values
The effectiveness of using optimised kernel parameters was tested on the b2 band and the low frequency band, which use the following set of dedispersion ranges:  

    b2 band:  
    range     0    150  0.1  1 1  
    range     150  300  0.2  1 1  
    range     300  500  0.25 1 1   
    range     500  900  0.4  2 2  
    range     900  1200 0.6  4 4  
    range     1200 1500 0.8  4 4  
    range     1500 2000 1.0  4 4  
    range     2000 3000 2.0  8 8  

    low frequency band:  
    range 0	  	 10.24	0.01 1 1  
    range 10.24		 50.56	0.01 2 2  
    range 50.56		 101.12	0.02 4 4  
    range 101.12	 202.24	0.04 8 8  

The optimal parameters found for each band are shown below:

![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/dm_range_parameters/params.png)

Optimised parameters for the B2 band. Black: Generic parameters, found by optimising over total speedup factor when all DM ranges are processed with the same parameters. Other: Optimised parameters for ranges 0-8.

![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/dm_range_parameters/low_params.png)

Optimised parameters for the low frequency band. Black: Generic parameters, found by optimising over total speedup factor when all DM ranges are processed with the same parameters. Other: Optimised parameters for ranges 0-3.

## Performance
The final speedup gained by using different optimised parameters for each DM ranges as opposed to the same parameters for all ranges was modest. As such, it is not recommended unless the code will be running unchanged on the same system for a long period of time. For the b2 band, the final speedup factor averaged over 5 trials was 6.28 with DM range specific optimisation vs 5.98 for generic parameters, an improvement in dpeedup factor of 5\%. This speedup factor includes time spent allocating global GPU memory and copying data to the GPU at the beginning of every time chunk. It was thought that the effect of individual kernel optimisation may be more pronounced for lower frequencies, where different DM ranges vary more significantly in aspect ratio. However, the results for this band were similar to the B2 band, with a speedup factor of 5.08 with DM range specific optimisation vs 4.75 for generic parameters, an again small improvement of 7.1\%.