# **Input File**

To run AstroAccelerate you must supply it with an input file that describes what you would like the code to do. This can contain a de-dispersion range, whether or not to mitigate for RFI, to zero-dm, to perform a single pulse search, output peaks in the processed data, perform a periodicity or Fourier domain acceleration search.  

To run the code with your input go to the "scripts" directory and then type:

_./astroaccelerate.sh ../input_files/your input file_

## **Example**

An example input file might look something like:

range	   0    150  0.1  1 1  
range	   150  300  0.2  1 1  
range	   300  500  0.25 1 1   
range	   500  900  0.4  2 2  
range	   900  1200 0.6  4 4  
range	   1200 1500 0.8  4 4  
range	   1500 2000 1.0  4 4  
range	   2000 3000 2.0  8 8  
sigma_cutoff	6  
sigma_constant  3.0  
max_boxcar_width_in_sec 0.5  
periodicity_sigma_cutoff 20  
periodicity_harmonics 32    
analysis  
-acceleration  
-output_ffdot_plan  
-output_fdas_list  
-periodicity  
-output_dmt  
-zero_dm  
-zero_dm_with_outliers  
-rfi  
-threshold  
-baselinenoise  
-fdas_custom_fft  
-fdas_inbin  
-fdas_norm  
debug  
-analysis_debug  
-failsafe  
file ~/filterbank/ska-mid-b2.fil  

The minus "-" in front of a keyword tells AstroAccelerate to ignore that input.

Lets look at the keywords and their arguments one by one.

_range_ tells the code to de-disperse and has input

range dm_low dm_high dm_step downsampling_in_input_time downsampling_in_output_time

so…

range 500  900  0.4  2 2

tells the code to de-disperse from a dm of 500 to a dm of 900 in steps of 0.4 (making (900-500)/0.4 dm trials) with input data downsampled (in time by 2, so 64uS would be binned into 128uS samples).

_sigma_cutoff_ is the SNR cutoff for your single pulse search. 

_sigma_constant_ is multiple of standard deviation which is used for outlier rejection

_max_boxcar_width_in_sec_ tells the single pulse search code what the maximum width (in seconds) of single pulses to search for.

_periodicity_sigma_cutoff_ is the SNR cutoff used in your search for periodic objects. Applies only when _periodicity_ is enabled.

_periodicity_harmonics_ in the number of harmonics used in harmonic summing.

_analysis_ tells the code to analyse the de-dispersed data outputting data into the output directory, these are binary files and so can be read in gnuplot (for example, you could use python etc) using:
splot "./57663_22588_B0531+21_000007.fil/global_analysed_frb.dat" binary format="%float%float%float%float" u 1:2:3 palette

_acceleration_ tells the code to perform a Fourier domain acceleration search on the data.  

_output_ffdot_plan_ tells the code to output the whole f f-dot plane from the FDAS search.  

_output_fdas_list_ tells the code to output just the peaks found in the f f-dot plane.  

_periodicity_ looks for periodic objects.    

_output_dmt_ outputs the entire de-dispersed data to a file (in ASCII).  

_zero_dm_ you can guess.  

_zero_dm_with_outliers_ this is part of the rfi mitigation routine and aims to remove RFI in spectra.

_rfi_ again this tries to eliminate RFI in channels.  

_threshold_ This outputs a thresholded de-dispersed plane (AstroAccelerate normally outputs peaks only).

_baselinenoise_ enables outlier rejection during calculation of mean and standard deviation.

_fdas_custom_fft_ Uses our fastest shared memory based FDAS search.  

_fdas_inbin_ Performs some interbinning out the output of FDAS.  

_fdas_norm_ De-reds the input to FDAS.  

_debug_ Outputs detailed information about the code and processing steps.  

_analysis_debug_ Uses a CPU code to perform analysis, this is used to debug our GPU code.

_failsafe_ Puts the code in failsafe mode. This turns off optimizations but ensures you get a result.

_file <your input file>_ This is the path to your input file.