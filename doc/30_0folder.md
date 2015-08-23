The 0/ folder
===============
The `0/folder` is required for every simulation and has to contain a certain amount of information for the initial conditions. One of the basic input from the user has to be `radius.json` where a spesific radius is assigned to every particle. 

Syntax
-------------------

````
    {
        "name": "radius",
        "data" :
        {
            "[particleID]": [ [number]],
            "[particleID]": [ [number]]
        }
    }
 
````

Example lines from `0/radius.json`
-----------------------------

````
    {
        "name": "radius",
        "data" :
        {
            "1": [ 3.0e-1],
            "2": [ 6.0e-1]
        }
    }
```` 

 Explanation
----------------
 - [particleID] is the number of particle (starting with `0`) the user defined
 - [number] is the wanted radius for a specific [particleID]
 

The following input is depending on the `<ModelEqnType>` and `<ModelEqnName>` defined along with the `modelEqn` command (compare `14_ModelEquation.md`). For every model Equation a file has to be created by the user which defines the initial profile. So for every grid point a specific value has to be assigned (compare `11_Mesh.md`). For 2 particles with a number of grid points equal to 6 the input of `<ModelEqnName>.json` would look like:

Syntax
-------------------

````
{
    "name": "`<ModelEqnName>`",
    "data" :
    {
        "1": [ [number], [number], [number], [number], [number], [number],[BCvalue]],
        "2": [ [number], [number], [number], [number], [number], [number],[BCvalue]]
    }
}
````

Example lines from e.g. `0/heat.json`
-----------------------------

````
{
    "name": "heat",
    "data" :
    {
        "1": [ 313.15, 313.15, 313.15, 313.15, 313.15, 313.15, 500],
        "2": [ 513.15, 513.15, 513.15, 513.15, 513.15, 513.15, 600]
    }
}
````
 Explanation
----------------
 - `<ModelEqnName>` is the name of the model Equation specified in in `in.file` (compare `14_ModelEquation.md`)
 - [number] is the wanted value for a specific [particleID]
 - [BCvalue] is the boundary value of the particle. The context of the value is defined after the boundary condition specified from the user for that specific model Equation (compare `14_ModelEquation.md`):
    * NEUMANN BC - the constant heat flux to the environment (could be changed via coupling model json!)
    * DIRICHLET BC - [BCvalue] not needed since the particle surface temperature stays constant
    * CONVECTIVE BC - the constant bulk property - e.g. the fluid temperature (could be changed via coupling model json!)
    * For transient boundary condition see `13_CouplingModel_JSON.md`
    
Read next
-----------
 - [40_testHarnessConfigFiles.md](40_testHarnessConfigFiles.md)

Related
----------
 - [11_Mesh.md](11_Mesh.md)
 - [13_CouplingModel_JSON.md](13_CouplingModel_JSON.md)
 - [16_EquationProperties.md](16_EquationProperties.md)
 - [30_0folder_1_PhaseFraction.md](30_0folder_1_PhaseFraction.md)

