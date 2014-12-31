Model Equation 1D Spherical
======================
This model equation allows the user to perform simulations using a one-dimensional discretization of a particle. Discretization is only in the radial coordinate of a spherical coordinate system located at the center of the particle.

The output of this model is a radial profile of, e.g., temperature or species concentration. The last (i.e., outermost) grid point is the temperature/concentration in the surrounding fluid.


Example lines from `in.file`
-----------------------------

````
   #Model Equations
   modelEqn 1DSpherical  <ModelEqnName>     BC0 1   BC1 2
```` 
 Explanation
----------------
 - `modelEqn` indicates that a modelEqn option is following
 - `1DSpherical` selects this model 
 - `<ModelEqnName>` defines the name of the equation. Currently implemented are `heat` and `species` are implemented.
 - To solve different e.g. `heat` equations just name them e.g. `heat1` and `heat2` 
 - `BC0` and `BC1` allow the user to define oundary conditions at the surface (`BC1`)  and center (`BC0`) of the particle. Note that boundary condition at the middle (`BC0`) will be overwritten by the assumption of symmetry! The user must specify a [number] which indicates the type:
    * 0 -> NEUMANN (fixed flux) **constant heat flux is specified as the last value in the initial conditions** (see 30_0folder.md)
    * 1 -> DIRICHLET (fixed value) 
    * 2 -> CONVECTIVE (convective transfer at the surface of the particle) **surrounding (fluid) temperature is specified as the last value in the initial conditions** (see 30_0folder.md)
 - See docu in `doc/pdf/1_modeleqn/ for more infos and underlying equations
   
 
Read next
-----------
 - [15_Chemistry.md](15_Chemistry.md)

Related
----------
 - [14_ModelEquation.md](14_ModelEquation.md)
