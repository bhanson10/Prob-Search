# Prob-Search
The "Prob-Search" repsository containes the codebase accompying the paper "An extensible framework for Probabilistic Search for stationary 
or randomly-moving targets characterized by super-Gaussian or experimentally-defined Regions of Interest (ROI)" published in ___. Below is an in-depth summary
explaining the proper usage of all the software included in said repository. In all, the codebase provides the computational framework necessary
to: <br> <br>
&ensp; &ensp; &ensp; 1. Generate the relaxation advection field that drives an arbitrary PDF to a steady-state shape for either analytical or numerical PDFs <br>
&ensp; &ensp; &ensp; 2. Drive said arbitrary PDF using the relaxation advection field and isotropic diffusion (as well as create animations of this process) <br>
&ensp; &ensp; &ensp; 3. Use the relaxaction advection field to perform a non-evasive search process with _n_ search vehicles on varying orbits <br>
&ensp; &ensp; &ensp; 4. Use the relaxaction advection field plus an additional field function to perform an evasive search process with _n_ search vehicles on varying orbits <br>

## relaxation_advection.m
The purpose of _relaxation_advection.m_ is to achieve Objective 1 and 2, generating the relaxation advection field for analytical and numerical PDFs and demonstrating that the relaxation advection field does as it should by driving an arbitrary PDF with it. We begin by generating an initial PDF to use as our test case. The following are possible options for initial PDFs: <br>
&ensp; &ensp; &ensp; 1. **Gaussian/Super-Gaussians** <br>
&ensp; &ensp; &ensp; &ensp; - Parameters that may be varied: <br>
&ensp; &ensp; &ensp; &ensp; &ensp; * _xvbar_: Mean of Gaussian/Super-Gaussian <br>
&ensp; &ensp; &ensp; &ensp; &ensp; * _P_: covariance of Gaussian/Super-Gaussian <br>
&ensp; &ensp; &ensp; 2. **Kidney Bean** <br>
&ensp; &ensp; &ensp; &ensp; - This PDF is a fairly rigid example, meant to demonstrate a PDF that is shaped like a kidney bean. For this reason, changing the default parameters may cause problems <br>
&ensp; &ensp; &ensp; 3. **Numerical PDFs** <br>
&ensp; &ensp; &ensp; &ensp; - This implementation takes a set of data points and creates a PDF out of them via the alpha-convex hull implementation or kernel density estimation. Change the _flag_ parameter to choose between the three provided datasets (information explaining each of the datasets throughouly can be found in the paper) or generate your own dataset via _create_dataset.m_. 
