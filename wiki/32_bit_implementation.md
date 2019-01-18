# **32-bit Implementation**

The code now reads 32 bit filterbank files. Because a significant re-write of the code would be needed to implement true 32 bit functionality we have produced a routine to scale the data to 16 bits. 
The data in the files are re-scaled to an unsigned short (16 bit value between 0 and 65535). This is achieved by...

1. Finding the maximum and minimum values in the input 32 bit filterbank files.
2. Next a binning value is calculated according to bin = (max - min) / 65535
3. Filterbank data are then read from the file and stored as an unsigned short using the following conversion:
16 bit value = (32 bit value - minimum) / bin

## **Results**

The following images show a comparison of two filterbank files created using sigproc fake. 
The only difference is one is 32 bit, the other 16 bit.

![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/32bit/one.png)
![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/32bit/two.png)
![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/32bit/three.png)
![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/32bit/four.png)
![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/32bit/five.png)

It's likely that the 32 bit values have slightly higher SNR recovery (approximately 10%) because they are re-scaled to the full range of the unsigned short (0 - 65535) and the 16 bit interbank is not using the full range. But this hasn't been verified.
For data that has true 32 bit fidelity we would expect a decrease in recovered SNR using this method as opposed to processing the data in its native 32 bits.