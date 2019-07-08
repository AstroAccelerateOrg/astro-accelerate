# **Fourier Domain Acceleration Search**

# Tests included in the AstroAccelerate FDAS code

There are two tests included with FDAS code. First test verifies the convolution code and second test aims to check whole FDAS algorithm. To perform these test ‘acceleration’ flag must be switched on in input text file and appropriate parameter must be #defined in `include/aa_params.hpp`. Some basic parameters of these test could be changed in `include/aa_fdas_test_parameters.hpp`. 

## Convolution code test
The convolution test consists from a signal which replaces data otherwise produced by the de-dispersion algorithm with a simple test signal which consists from number of top-hat signals of different width and saw-tooth wave. This signal is written into a file 'acc_conv_test_input_signal.dat'. The aim of this test is to verify that convolution result produced by our implementation of overlap-and-save is correctly producing continuous linear convolution. This test signal is then convolved with series of top-hat filters of increasing width. Result of this test is then exported into the file ‘acc_fdas_conv_test.dat’. After exporting this file, code will exit. No data from filterbank file are processed. To enable this test a parameter `FDAS_CONV_TEST` must be #defined in ‘include/aa_param.hpp’. Parameters of this test are `FDAS_TEST_FILTER_INCREMENT` which controls how width of the top-hat filters grow with given template index and `FDAS_TEST_TOOTH_LENGTH` which controls width of tooths in saw-tooth wave (both parameters defined in ‘include/aa_fdas_test_parameters.hpp’). Output of this test is shown in the figure below.

![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/FDAS/simple_convolution_test.png)
![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/FDAS/tophat_signal.png)

This could be plotted using gnuplot:
* `set xrange[49.5:50.5]`
* `set yrange[0:30000]`
* `splot 'acc_fdas_conv_test.dat' using 1:2:3 palette`


## Accelerated sinusoid signal
To test the whole FDAS algorithm we generate accelerated sinusoid signal according to

![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/FDAS/f1.png)

![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/FDAS/f2.png)

![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/FDAS/f3.png)

![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/FDAS/f4.png)


Where *t* is time, *a* is pulsar acceleration, *a* is speed of light, *T* is observation length, *f* is pulsar period, *d* is pulsar duty cycle, *H* is number harmonics, *A* is amplitude and *Z* is acceleration template. Parameters which could be changed in the ‘include/aa_fdas_test_parameters.hpp’ are listed together with default value here:

* *f*: `FDAS_TEST_FREQUENCY 105.0`
* *Z*: `FDAS_TEST_ZVALUE 6`
* *H*: `FDAS_TEST_HAMONICS 4`
* *d*: `FDAS_TEST_DUTY_CYCLE 1.0`
* *A*: `FDAS_TEST_SIGNAL_AMPLITUDE 1.0`

Output of this test has the same format as normal FDAS algorithm output. Again code will exit after exporting result of the test into the file. To enable this test a parameter `FDAS_ACC_SIG_TEST` must be #defined in ‘include/aa_param.hpp’. Output of this test is shown in the figure below.

![](https://github.com/AstroAccelerateOrg/images/blob/master/wiki/FDAS/acceleration_test_plane.png)
