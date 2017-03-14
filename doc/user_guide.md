== User Guide

= Basic

This code should work on any Fermi, Kepler or Maxwell GPU. 
CUDA version 8.0 and above is needed. 

This code uses one GPU. 
If you have a multi-gpu system, use the command "nvidia-smi" to figure out which card you want to use.

$ nvidia-smi
Mon Dec  5 10:08:56 2016       
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 367.44                 Driver Version: 367.44                    |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|===============================+======================+======================|
|   0  GeForce GTX 1080    Off  | 0000:01:00.0     Off |                  N/A |
| 27%   39C    P8     6W / 180W |      1MiB /  8113MiB |      0%      Default |
+-------------------------------+----------------------+----------------------+
|   1  GeForce GTX 980     Off  | 0000:06:00.0      On |                  N/A |
|  0%   59C    P0    50W / 185W |    358MiB /  4028MiB |      0%      Default |
+-------------------------------+----------------------+----------------------+
                                                                               
+-----------------------------------------------------------------------------+
| Processes:                                                       GPU Memory |
|  GPU       PID  Type  Process name                               Usage      |
|=============================================================================|
|    1      6789    G   /usr/lib/xorg/Xorg                             269MiB |
|    1      7317    G   compiz                                          86MiB |
+-----------------------------------------------------------------------------+

If you want to use the GTX 1080, simply edit the file "params.h" located in lib/AstroAccelerate and set the CARD variable to 0:
#define CARD 0

= Compilation
- Using Makefile:
Go to the "lib" directory, and use the command "make" to create the executable:
$ make clean && make -j 8
(this assumes that cuda and drivers are installed)
- Using CMake:
$ cd build
$ cmake .. # This will generate a Makefile
$ make 	 # make the code and create the executable

= Test files

Create test files, such as:

fake -nbits 32 -tobs 60 -dm 45 -tsamp 491.52 -fch1 151 -foff 0.012207031250 -nchans 2592 -period 500 -snrpeak 0.5 > lotaas_500_0.5.dat
fake -nbits 32 -tobs 60 -dm 45 -tsamp 491.52 -fch1 151 -foff 0.012207031250 -nchans 2592 -period 500 -snrpeak 1 > lotaas_500_1.dat
fake -nbits 32 -tobs 60 -dm 45 -tsamp 491.52 -fch1 151 -foff 0.012207031250 -nchans 2592 -period 500 -snrpeak 8 > lotaas_500_8.dat
fake -nbits 32 -tobs 60 -dm 45 -tsamp 491.52 -fch1 151 -foff 0.012207031250 -nchans 2592 -period 18.125 -snrpeak 0.25 > lotaas_18.125_0.25.dat

= Input file

To run the code, you need an input file which looks like this: 

range   0    60    0.050 1  1
range   60   120   0.050 2  2
range   120  240   0.100 4  4
range   240  480   0.200 8  8
sigma_cutoff	7
analysis
acceleration
periodicity
output_dmt
zero_dm
rfi
debug
file <location of your filterbank file>

To switch a feature off, add a "-" in front of the associated keyword.
Here are some explanations about the keywords:

range <start dm> <end dm> <dm step size> <decimation in time, in f,t input> <decimation in time, in dm,t output>:
The above example "bins" at about 6x the diagonal dm in this example.

sigma_cutoff <value>:
This is the number of sigma (SNR) above which candidates will be reported. Currently this is the same for single pulses and periodicity search.

power <value>: 
this changes the 1/f^2 to 1/f^power in the dm trial.

analysis:
This turns on and off the single pulse search.

acceleration:
[Insert appropriate description]

periodicity:
This looks for periodic things.

output_dmt:
This outputs the de-dispersed space into a ASCII text file.

zero_dm:
[Insert appropriate description]

rfi:
[Insert appropriate description]

debug:
This outputs lots of info

file:
This is the location of the input file.

= Select the best parameters

A script "profiling.sh" is provided. It runs the code several times with different parameters, and select the optimal configuration. It then writes this configuration in the file "params.h".  

= Run the code

Go to the script directory and type the following command to run the code:
$ ./astro-accelerate.sh ./lotaas_18.txt

