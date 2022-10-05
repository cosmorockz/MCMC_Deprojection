# MCMC_Deprojection
X-ray spectrum analysis code written in python and C

The first part of our code is to create fake (.pha and .rsp files) from spectra provided in .dat/.spec/.txt files.
Then we divide our spectra into annuli and deproject them. We fit them using XSpec.
We then take the fitted spectrum from the projected and deprojected models and pass them to the next step.

In the next step, we take the density and temperature profiles and fit them using non-linear fitting
algorithms and MCMC to obtain G200, c and v_c parameters. This parameters then can be used to 
obtain the cooling time, free time and gravity profiles.
