# asteroseismic_fitting
Homemade tools to perform the matching of periods in asteroseismic fitting. Designed to work with pulsation output from the WDEC (see repository by that name)

The folder mkgrid contains simple code and a sample input file to generate lists of parameters to be ran through WDEC.
It generates the gridparameters file.

The folder fitper contains code to perform period matching. It is quite hands on to use and is for the advanced Fortran user.
It can serve as a starting point to any list matching you may be trying to do. It is specifically geared to be used with 
output from the WDEC. It requires two input files to work: calcperiods (output from WDEC) and obsperiods (input from the user).
A sample of each is provided if you want to run the code.

The code identifies the periods listed in calcperiods that best match the periods provided by the user in obsperiods.
As a start, look at two source files for information: the header of main.f90, and line 136-138 of fit.subroutine.f90.
The former explains (inperfectly) how to list periods in obsperiods while the latter shows how the fitness parameter (which
describes how good the match is) is computer. Weights can be assigned to each periods in obsperiods. In the sample file, all
periods are weighted equally (all have a weight of 1).

A handy feature for debugging purposes is the logical parameter called "verbose" in modules.f90. It shows details of how the fitting is done. It's good to look at that output on a single list of periods (the sample calcperiods only has one model and so is a great sample to use) and check that the fitting is done correctly before launching the code on a series of models.

The model parameters and corresponding matched periods are written out to the file fitnessparameters.dat.
