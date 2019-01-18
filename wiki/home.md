![](http://www.oerc.ox.ac.uk/sites/default/files/uploads/ProjectFiles/AstroAccelerate/AstroAccelerate.jpg)

# **Home**

## **Welcome to AstroAccelerate!**

## **Table of Contents**
* [Home](home.md)
* [32-bit Implementation](32_bit_implementation.md)
* [Calculation of mean and standard deviation](calculation_of_mean_and_standard_deviation.md)
* [Dedispersion kernel parameter optimisation](dedispersion_kernel_parameter_optimisation.md)
* [Fourier Domain Acceleration Search](fourier_domain_acceleration_search.md)
* [Past and present contributors to AstroAccelerate](past_and_present_contributors_to_astroaccelerate.md)
* [RFI mitigation](rfi_mitigation.md)
* [Single pulse search](single_pulse_search.md)
* [Using input files](using_input_files.md) 

## **Introduction**

AstroAccelerate is a many core accelerated software package for processing time domain radio-astronomy data. In our git repo you will find:  

1. A standalone code that can be used to process filterbank data.
2. A library that can be used to enable GPU accelerated single pulse processing (SPS) or Fourier Domain Acceleration Searching (FDAS).

### **Time domain radio-astronomy**

Pulsars are fast-spinning neutron stars which emit radio beams along their magnetic axis. Since they are rotating, the radio emissions - or pulses - are observed periodically from earth. On a large timescale these pulses are very precise, as such pulsars can be thought of as "cosmic clocks". The physics of pulsars is extreme: gravity, density and magnetic fields and so can't be reproduced in a laboratory. Therefore, study of pulsars can help physicists to gain new knowledge of the universe.

### **AstroAccelerate**

AstroAccelerate is a GPU-enabled software package that focuses on enabling real-time processing of time-domain radio-astronomy data. It uses the CUDA programming language for NVIDIA GPUs.  
The massive computational power of modern day GPUs allows the code to perform algorithms such as dedispersion, single pulse searching and Fourier Domain Acceleration Searching in real-time on very large data-sets which are comparable to those which will be produced by next generation radio-telescopes such as the SKA.  