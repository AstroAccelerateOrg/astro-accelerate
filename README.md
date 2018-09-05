Repository for Astro-Accelerate Code
====

Introduction
===
This documentation summarises the content of the Astro-Accelerate software.

Features
===
Astro-Accelerate is useful for a number of real-time astronomy applications.
Features include:
Acceleration Searching
Zero DM
RFI Mitigation

Software Inputs and Outputs
===
The software input is a sample data file.
To process the data file, astro-accelerate makes use of a configuration file.
Please see the section `Creating An Input Configuration File` for instructions on how to create a configuration file to process an input data file.
The software output is dependent on the choice of analysis module that is run.

Checking the Configuration of the Graphics Processing Unit (GPU) and Support for CUDA

In a terminal window, type

    nvidia-smi

The output will look similar to the following example
```
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 375.20                 Driver Version: 375.20                    |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|===============================+======================+======================|
|   0  Tesla P100-PCIE...  On   | 0000:02:00.0     Off |                    0 |
| N/A   32C    P0    26W / 250W |      0MiB / 16308MiB |      0%      Default |
+-------------------------------+----------------------+----------------------+
|   1  Tesla P100-PCIE...  On   | 0000:03:00.0     Off |                    0 |
| N/A   34C    P0    25W / 250W |      0MiB / 16308MiB |      0%      Default |
+-------------------------------+----------------------+----------------------+
                                                                              
+-----------------------------------------------------------------------------+
| Processes:                                                       GPU Memory |
|  GPU       PID  Type  Process name                               Usage      |
|=============================================================================|
|  No running processes found                                                 |
+-----------------------------------------------------------------------------+
```

Selecting the Graphics Processing Unit (GPU) in the case of a system with more than one GPU
===
If you have a multi-GPU system, you need to modify two files to use one specific card.

1. Edit `lib/AstroAccelerate/params.h`, so that the `CARD` variable is set to the right card number id:
    #define CARD ... 
2. Edit the Makefile and set the `GENCODE_FLAG` to the correct architecture. For example:
    GENCODE_FLAGS := $(GENCODE_SM61)
    Where `GENCODE_SM61` corresponds ot the Pascal architecture.

Software Pre-Requisites
===
CUDA: CUDA 8.0 (see https://developer.nvidia.com/cuda-downloads)

C++:

Compiler: (requires OpenMP support)

Creating An Input Configuration File
===
An example input configuration file is shown below:
 
    range    0    150  0.1  1 1
    range    150  300  0.2  1 1
    range    300  500  0.25 1 1
    range    500  900  0.4  2 2
    range    900  1200 0.6  4 4
    range    1200 1500 0.8  4 4
    range    1500 2000 1.0  4 4
    range    2000 3000 2.0  8 8
    sigma_cutoff     6
    analysis
    -acceleration
    -periodicity
    -output_dmt
    -zero_dm
    -zero_dm_with_outliers
    -rfi
    fdas_custom_fft
    -fdas_inbin
    fdas_norm
    debug
    file /home/wa78/filterbank/ska-mid-b2.fil`
 
Features can be turned on or off by adding a character at the beginning of the line (here "-" is used).

Setting `range`
=
`range` tells the code to dedisperse and has input.

The format for the `range` parameter is
`range dm_low dm_high dm_step downsampling_in_input_time downsampling_in_output_time`
where `dm_low`, `dm_high`, `dm_step`, `downsampling_in_input_time`, and `downsampling_in_output_time` are to be replaced with suitable numerical values.
For example, a valid input for the `range` parameter is

    range 500  900  0.4  2 2
 
astro-accelerate will parse this input and dedisperse from a `dm` of `500` to a `dm` of `900` in steps of `0.4` (making (`900`-`500`)/`0.4` `dm` trials), with input data 
downsampled in time by `2` (e.g. 64 uS would be binned into 128 uS samples).
 
Setting `sigma_cutoff`
=
`sigma_cutoff` is the Signal-to-Noise Ratio (SNR) cutoff for your single pulse search.

Setting `analysis`
= 
`analysis` this tells the code to analyse the dedispersed data, outputting data into the output directory. The output data are binary files and so which can be read in gnuplot and python.

For example, you could use python as follows
    splot "./57663_22588_B0531+21_000007.fil/global_analysed_frb.dat" binary format="%float%float%float%float" u 1:2:3 palette
 
 Setting `acceleration`
 =
`acceleration` this does a Fourier domain acceleration search on the data.

Setting `periodicity`
=
`periodicity` sets a search for periodic objects.

Setting `output_dmt`
=
`output_dmt` outputs the entire dedispersed data to a file (in ASCII).

Setting `zero_dm`
= 
`zero_dm` you can guess (well, I can't!)

Setting `zero_dm_with_outliers`
=
`zero_dm_with_outliers` is part of the RFI mitigation routine.
 
Setting `rfi`
=
`rfi` tries to eliminate RFI.
(Astro-Accelerate welcomes developers to supply the team with data that includes RFI, which would be very helpful.)

Setting `fdas_custom_fft`
=
`fdas_custom_fft` runs fourier domain acceleration search with a custom FFT.

Setting `fdas_inbin`
=
`fdas_inbin` performs interbinning on the complex output.

Setting `fdas_norm`
=
`fdas_norm` performs PRESTO block median normalization.

Setting `debug`
`debug` this gives detailed output.
 
Setting `file`
Please supply the input data file by using the `file` setting followed by the path to a valid
input data file

    file <your input file>
 
Input data files are by definition assumed to be formatted with 8-bit data.
Please contact Astro-Accelerate if you have input data that is not 8-bits let me know (development for support for this is in progress and we may be able to help).

Installation and Usage (Linux)
===
Step 1: Download the files
==
Obtain the Astro-Accelerate code by doing
    git clone https://github.com/AstroAccelerateOrg/astro-accelerate

Step 2: Compile
==
Ensure you have the correct environment and pre-requisites.
Set-up the environment

    ./bootstrap.sh

Next, `cd` into the `lib` directory by doing

    cd lib/

The software can be compiled using the supplied makefile. To do so, simply type

    make

The compilation process indicates which components are being compiled. There should be no errors or warnings.

The result will be an executable called

    astro-accelerate

in the astro-accelerate directory.

Step 3: Run
==
Astro-Accelerate assumes its input is ready and compatible. To obtain compatible input, please follow the steps below. Please also ensure that the aforementioned environment variables have been set.

1. Run astro-accelerate using the format

    ./astro-accelerate --settingname settingvalue

By default, the output of astro-accelerate will be located in the same directory in which astro-accelerate was executed.
Configuration files may be used to further specify, set, and change options.
By default, the astro-accelerate executable looks for a configuration file in the same directory as the executable.
A configuration file is required in order to run astro-accelerate.
A number of example configuration files are included in the repository.

Additional Options
==
All settings in the configuration file can be overriden in the CLI using:

    --settingname value

where `settingname` is the exact name of the setting (case-sensitive) in the configuration file, and `value` is a suitable valid value of the parameter.

Failure to provide a valid input will stop execution. It is recommended to leave all settings
as specified by default.

A custom configuration path can be specified using the following option

    --astro-accelerate-config-path /path/to/config.txt

Optimisations and Tuning
==
Astro-Accelerate comes with the facility to tune the software to the input that the user provides.
1. To do this, `cd` to the `scripts` directory.
Modify `profiling.sh`, changing the line that says

        ./astro-accelerate.sh ../input_files/ska_tune.txt

    to

        ./astro-accelerate.sh ../input_files/<your input file>

    where `<your input file>` must be replaced by the input configuration file.

3. The next step is to run the profiler tool that is provided in the repository

        ./profiling.sh

    This will create an optimised code for your search and GPU type.

4. Then, astro-accelerate can be run as usual by doing 

        ./astro-accelerate.sh ../input_files/<your input file>

    where `<your input file>` must be replaced by the path to the input configuration file as specified in the previous step.

Further Documentation
===
More detailed information can be found on the Wiki page of the repository
https://github.com/AstroAccelerateOrg/astro-accelerate/wiki

and the Astro-Accelerate webpage
http://www.oerc.ox.ac.uk/projects/astroaccelerate

Contact and Support
===
If you notice any errors or have suggestions,
please contact someone on the astro-accelerate team,
or file an issue or bug report on the repository.

Disclaimers and Licensing
===
A number of code files may be covered by different licences.
Please refer to these seperately.
Please also refer to the Astro-Accelerate licence file provided.

Copyright © 2018 Astro-Accelerate. All rights reserved.
