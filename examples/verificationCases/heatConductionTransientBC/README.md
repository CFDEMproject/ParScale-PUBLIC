Heat Conduction with Transient Boundary Conditions
=======================================================

This a checked test case for the transient change of boundary conditon in comparism with an analytical solution.

How to run:
-----------

- Open a terminal at the `in.file_*`´s postion and just type `./Allrun_convective`
- In case `Allrun_*` and `Allclean` are not executable type `chmod +x All*`
- After the simulation is done type `./Allclean` to delete all generated files

Description:
------------

- Cooling of a sphere (with temperature specified in '0/heat.json') in a fluid with a temperature ramp specified in the 'settings/coupling_json.json' file.
- A plot is prepared that shows the ambient (fluid) temperature, as well as that of the particle at certain locations (see the 'plotMe.m' file)
