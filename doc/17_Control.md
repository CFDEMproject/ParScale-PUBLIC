Control Command
======================

The Control Command let the user decide how long a simulation will run. Also it gives the options to adjust time step and output time step.


Syntax
-------------------

````
    #Control
    control outputTimeStep [number]
    control timeStep [number]
    control run [number]
````

Example lines from `in.file`
-----------------------------

````
    #Control
    control outputTimeStep 0.02
    control timeStep 0.01
    control run 0.050
```` 

 Explanation
----------------
- `outputTimeStep` defines the interval after each data output is written
- `timeStep` defines the absolute timestep for the integrator (automatic sub time stepping)
- `run` defines the absolute time the simulation will run
- Note that `timeStep` and `run` have to be deleted in case of a coupling model to LIGGGHTS which runs as the master program. 

Read next
-----------
 - [18_PhaseChange.md](18_PhaseChange.md)


Related
----------
 - [13_CouplingModel_LIGGGHTS.md](13_CouplingModel_LIGGGHTS.md)


