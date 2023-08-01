# Prob-Search
The "Prob-Search" repsository containes the codebase accompying the paper "An extensible framework for Probabilistic Search for stationary 
or randomly-moving targets characterized by super-Gaussian or experimentally-defined Regions of Interest (ROI)" published in ___. Below is an in-depth summary
explaining the proper usage of all the software included in said repository. In all, the codebase provides the computational framework necessary
to: <br> <br>
&ensp; &ensp; &ensp; 1. Generate the relaxation advection field that drives an arbitrary PDF to a steady-state shape for either analytical or numerical PDFs <br>
&ensp; &ensp; &ensp; 2. Drive said arbitrary PDF using the relaxation advection field and isotropic diffusion (as well as create animations of this process) <br>
&ensp; &ensp; &ensp; 3. Use the relaxaction advection field to perform a non-evasive search process with _n_ search vehicles on varying orbits <br>
&ensp; &ensp; &ensp; 4. Use the relaxaction advection field plus an additional field function to perform an evasive search process with _n_ search vehicles on &ensp; &ensp; &ensp; &ensp; varying orbits <br>

## relaxation_advection.m
The purpose of _relaxation_advection.m_ is to achieve Objective 1 and 2, generating the relaxation advection field for analytical and numerical PDFs and demonstrating that the relaxation advection field does as it should by driving an arbitrary PDF with it. 
