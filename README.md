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

&ensp; &ensp; &ensp; 1. **Gaussian/super-Gaussians** <br>
&ensp; &ensp; &ensp; &ensp; - Parameters that may be varied: <br>
&ensp; &ensp; &ensp; &ensp; &ensp; * _xvbar_: Mean of Gaussian/super-Gaussian <br>
&ensp; &ensp; &ensp; &ensp; &ensp; * _P_: covariance of Gaussian/super-Gaussian <br>
&ensp; &ensp; &ensp; &ensp; &ensp; * _s_: order of Gaussian/super-Gaussian ($s>1$ is considered a super-Gaussian in this context)<br>
&ensp; &ensp; &ensp; 2. **Kidney Bean** <br>
&ensp; &ensp; &ensp; &ensp; - This PDF is a fairly rigid example, meant to demonstrate a PDF that is shaped like a kidney bean. For this reason, changing the default parameters may cause problems <br>
&ensp; &ensp; &ensp; 3. **Numerical PDFs** <br>
&ensp; &ensp; &ensp; &ensp; - This implementation takes a set of data points and creates a PDF out of them via the alpha-convex hull implementation or kernel density estimation. Change the $flag$ parameter to choose between the three provided datasets (information explaining each of the datasets throughouly can be found in the paper) or generate your own dataset via _create_dataset.m_. <br>

Some other paramters that are available for changing are: <br>
&ensp; &ensp; - _dt_: time step of the simulation <br>
&ensp; &ensp; - _N_: number of grid cells in one dimension of the $N\times N$ grid <br>
&ensp; &ensp; - _L_: length and width of the grid <br>
&ensp; &ensp; - _lambda_: homogenous isotropic diffusion constant, such that $D=\lambda I$ <br><br>

To choose one of these specific PDF options, comment out the other instances that call new PDFs. Once the initial PDF is chosen, the code generates the relaxation advection field in the $x$ and $y$ directions, each an $N\times N$ grid with the advection value at the grid center. The next step generates a uniform PDF of the same size _PDF_U_ and propagates that PDF using the relaxation advection field and homogeneous isotropic diffusion until convergence (convergence is explain in the paper). The _create_video_ function combines all of the frames generated by the code into a single .MP4.

## non-evasive_search.m

The non-evasive search code demonstrates how the relaxation advection can be used to demonstrate the rate of change of confidence in observations searching for a stationary target in a ROI. For the non-evasive target, we assume that the location of the search vehicles has no impact on the behavior of the target, thus the observations act as negative driving forces that lower the probability at the location of the search vehicle, and the relaxation advection acts as a driving force towards the steady-state. After the observations are done being taken, as $t\rightarrow \infty$ the PDF returns to its steady-state. 

### Tunable Simulation Constansts
&ensp; &ensp; - _.T_: total time of the search <br>
&ensp; &ensp; - _.dt_: time step size of simulation <br>
&ensp; &ensp; - _N_: size of cartesian grid <br>
&ensp; &ensp; - _L_x_,_L_y_: length of cartesian grid in $x$-, $y$-direction <br>

### Tunable Target Constansts
&ensp; &ensp; - _.lambda_: same as above <br>
&ensp; &ensp; - _.stats_flag_: choose which initial PDF to start with <br>
&ensp; &ensp; &ensp; &ensp; * _stats_flag_$==1$: Gaussian statstics described by parameters above <br>
&ensp; &ensp; &ensp; &ensp; * else: numerical statistics dependent on _flag_ where the datasets are the same as previously stated <br>

### Tunable Drone/Search vehicle Constansts
&ensp; &ensp; - _.num_: number of search vehicles <br>
&ensp; &ensp; - _.ang_speed_: andgular speed of each search vehicle <br>
&ensp; &ensp; - _.init_theta_: initial starting position of each search vehicle <br>
&ensp; &ensp; - _.d_: disruptivity of each drone (see paper) <br>
&ensp; &ensp; - _.orbit_flag_: flag that dictates which orbit path is taken <br>
&ensp; &ensp; &ensp; &ensp; * _orbit_flag_$==1$: Concentric circles where radius can be varied <br>
&ensp; &ensp; &ensp; &ensp; * _orbit_flag_$==2$: Cassini ovals with tunable parameters <br>
&ensp; &ensp; &ensp; &ensp; * _orbit_flag_$==3$: Lemniscates with tunable parameters <br>
&ensp; &ensp; - _.sigma_: width of field-of-view of observations of search vehicle <br>

The simulation generates the initial PDF chosen as well as the orbits of the search vehicles, then demonstrates how the PDF changes as the observations drive the probabilty down and the advection field drives the proabability back to its steady-state. Again, _create_video_ can be used to save the frames in a .MP4 file. 

## evasive_search.m

The evasive search code demonstrates how the relaxation advection can be used to demonstrate the rate of change of confidence in observations searching for a stationary target in a ROI. For the evasive target, we assume that the magnitude and direction of the target is impacted by the locations of the search vehicles, and that the target's velocity will point away from the center of the search vehicle. If the search vehicle is closer by, the magnitude of the velocity will be higher, and the random motion will be more sporadic. The evasive velocity is modeled as the negative gradient of a super-Gaussian with $s=0.6$ (so as to avoid increasingly larger velocities if the search vehicle and the target share the same location) and the diffusion is modeled as the superposition of super-Gaussians placed at the locations of the search vehicles. In the background of the total advection field is the relaxation advection, such that there is still a desire to return to the general ROI, though it is scaled so as to not be overbearing. 

### Tunable Target Constansts (not already stated)
&ensp; &ensp; - _.psi_: skittishness of target (see paper) <br>


### Tunable Drone/Search vehicle Constansts (not already stated)
&ensp; &ensp; &ensp; &ensp; * _orbit_flag_$==4$: Rotating, "herding" Cassini ovals, preferable orbit family for evasive search (see paper) <br>
&ensp; &ensp; - _.sigma_: width of field-of-view of observations of search vehicle <br>
&ensp; &ensp; - _.P_: covariance of super-Gaussian used for determining advection and diffusion <br>

The simulation generates the initial PDF chosen as well as the orbits of the search vehicles, then demonstrates how the PDF changes as the observations drive the probabilty down and the advection field drives the proabability around depending on the location of the drones as well as the steady-state. An apparent herding effect is demonstrated when using the "herding" orbit path, as is the goal for evasive target searching. <br> <br>

For further information about code usage, please contact blhanson@ucsd.edu
