Verification case one-dimensional Stefan-Diffusion
=======================================================
This is a test case to verify the correct implementation of a three-phase porous particle including evaporation based on wet core evaporation. The analytical solution is presented in `Analytical_solution.pdf`. For the calculation see `Main_stefan_diffu.m`,`evaporation_time.m`, `StefanDiffusion1DSpherical.m` and `diameter_evo.m` in the `octave/` folder.

How to run:
-----------
- Open a terminal at the `in.file_*`´s postion and just type `./Allrun_convective`. 
- In case `Allrun_*` and `Allclean` are not executable type `chmod +x All*`
- After the simulation is done type `./Allclean` to delete all generated files

Description:
------------
- Evaporation of wet core in the particle
- Isothermal conditions
- Constant solid phase fraction
- Ethanol as liquid phase 
- `Allrun` will execute 'ParScale' and the analytical solution
- Plotting functionality for different times is provided in `octave/plot_liquid.m` and can be changed easily

