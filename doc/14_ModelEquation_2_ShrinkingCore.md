Model Equation Shrinking Core
======================
This model equation allows the user to perform simulations using a simple shrinking core model. 

The primary output of this model is the position of the reaction front (made dimensionless with the particle radius). The second grid point is the temperature/concentration in the surrounding fluid.


Example lines from `in.file`
-----------------------------

````
   #Model Equations
   modelEqn ShrinkingCore  <ModelEqnName>    BC1 2
```` 
 Explanation
----------------
 - `modelEqn` indicates that a modelEqn option is following
 - `ShrinkingCore` selects this model 
 - `<ModelEqnName>` defines the name of the equation. Currently, the shrinking core model can be only used to study (reactive) mass transfer, hence the name must be `species`.
 - `BC1` allows the user to define boundary conditions at the surface (`BC1`) of the particle. The user must specify a [number] which indicates the type:
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
