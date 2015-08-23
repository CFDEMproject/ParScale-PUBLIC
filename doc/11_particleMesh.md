Particle Mesh 
======================
The mesh information given by the user defines the number of grid points for which one particle is discretized.

Syntax
-------------------
````
    #Mesh
    particle_mesh nGridPoints [number]
````

Example lines from `in.file`
````
    #Mesh
    particle_mesh nGridPoints 9 
````


Explanation
------------
- `particle_mesh nGridPoints` equals the number of grid points

The radius for every particle must be specified in the file `0/radius.json` (units in [m])

Syntax
-------------------
````
    {
        "name": "radius",
        "data" :
        {
            "1": [ [number] ],
            "2": [ [number] ]
        }
    }
````

Example lines from `0/radius.json`
-------------------

````
    {
        "name": "radius",
        "data" :
        {
            "1": [ [5e-3] ],
            "2": [ [5e-4] ]
        }
    }
````

Explanation
------------
- Radius for every particle can be set according to the Particle ID (0,1,2,..) 
- Note that in case of a coupling to LIGGGHTS this input is overwritten by the pulled values from LIGGGHTS.
 

Read next
-----------
 - [12_particleData.md](12_particleData.md)

Related
-----------
 - [13_CouplingModel.md](13_CouplingModel.md) (for more information on coupling models)