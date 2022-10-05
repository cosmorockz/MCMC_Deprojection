# date: 05.11.2019

# This folder contains the following
1. avg_spectra.c : to calculate the averaged spectra by binning the projected
    spectra from pass_1.3/. The inputs are redshift (1e-3, default) and the 
    projected spectra in pass_1.3/apec_.../. 
    Run it with (as and example)
    ./avg_spectra 0.001 ../../pass_1.3/apec_0100/ 
     
2. In the next step, we want to fit the binned spectra and deproject them to
   obtain the density, temperature profiles. Ask Arjun to do this. 

3. We then want to combine the individual density and temperature profiles 
   into a single file to feed into the fitting machine. 

4. The fitting machine is compiled using 'make' 

5. it is then run with
   ./main run_nulsen_fit:yes/no 
   Remember to set 'run_nulsen_fit:no' if you do not need run the MCMC fit 
   of nulsen. It takes quite some time to finish the fit. 

6. Completed. 
