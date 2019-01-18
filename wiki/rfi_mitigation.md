# **Basic RFI Mitigation**

AstroAccelerate has some basic routines that can be used to remove RFI.  
There is a standard zero-dm option which is enabled using the keyword "zero_dm" in the input file.  
There is an option to perform zero-dm with outlier rejection "zero_dm_with_outliers".  
There is an option to perform outlier rejection on individual channels "rfi".  

## **Results**

The following gif shows images of data collected on the Lovell telescope<sup>1</sup> processed using AstroAccelerate.  

Image Left:   No RFI mitigation.  
Image Center: Old RFI AstroAccelerate Algorithms.  
Image Right:  New algorithms using a moving average (enabled with both "zero_dm_with_outliers" and "rfi" keywords).  

![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/rfi_removal/out_ani.gif)

<sup>1</sup> With great thanks to Mitch Mickaliger and Ben Stappers for finding test data for us to use.  