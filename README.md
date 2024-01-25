# Deconvolution
Matlab codes to deconvolve two 1-D signals (i.e., time series)

As written, the matlab reads in a .mat file that has 3 column vectors of equal length: time, input, and output.  There are example .mat files included here.

The code is based on the algorithm of Cirpka et al., (2007).  The original code was written by Olaf Cirpka.  I added a routine to discern and use the sample autocovariance function of the transfer function in order to estimate the transfer function itself, on the next iteration.  

The data included in the several examples here are from a tracer test described in Gambel et al., Water Resources Research, in press.  Please contact me at dbenson@mines.edu for an update on that citation.

In the folders you'll find the .m files "deconv_dave_2.m" that look for a specific .mat file in the folde and perform the deconvolution. The user must specify the name of the .mat file with input data, the distance between input and output signals (assumes collection in a stream, say), and an initial guess at the covariance function of the filter, and a maximum amount of epistemic noise.  There is also a data_prep.m file that you can use to make your .mat files.  
