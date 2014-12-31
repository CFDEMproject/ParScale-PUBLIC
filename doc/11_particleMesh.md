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
            "0": [ [number] ],
            "1": [ [number] ]
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
            "0": [ [5e-3] ],
            "1": [ [5e-4] ]
        }
    }
````

Explanation
------------
- Radius for every particle can be set according to the Particle ID (0,1,2,..) 
- Note that in case of a coupling to LIGGGHTS this input is overwritten by the pulled values from LIGGGHTS.
 

Read next
-----------
 - [12_ParticleNumber.md](12_ParticleNumber.md)

Related
-----------
 - [13_Coupling.md](13_Coupling.md) (for more information on coupling models)