Heat Conduction with Transient Boundary Conditions
=======================================================

This a checked test case for the transient change of boundary conditon in comparism with an analytical solution. The case considers unsteady heat conduction in a multi-component sphere (made up of solid, liquid, and gas).

How to run:
-----------
- Make all settings in the "octave/preProcessor.m" file. In case the user would like to compare to an analytical solution, the 'evaporationRateConstant' must be set to zero. Otherwise, evaporation will take place, not allowing for an analytical solution.
- Open a terminal at the `in.file_*`´s postion and just type `./Allrun_convective`. This will execute the preProcessor.m file (and hence make settings in the JSON files), and then start the ParScale simulation
- In case `Allrun_*` and `Allclean` are not executable type `chmod +x All*`
- After the simulation is done type `./Allclean` to delete all generated files

Description:
------------
- Convective heating/cooling of a sphere with no coupling (default) in a surrounding gas with temperature as specified in '0/heat.json'.
- In case JSON-coupling is used: heating/cooling of a sphere with a temperature ramp specified in the 'settings/coupling_json.json' file.
- A plot is prepared that shows the ambient (fluid) temperature, as well as that of the particle at certain locations (see the 'plotMe.m' file). This is deactivated by default.
