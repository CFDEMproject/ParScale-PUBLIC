Particle Data
======================
The particle number information given by the user defines the number of particles that are considered in one simulation.

Syntax
-------------------
````
#Particle Number
particle_data number_particles [number]
````

Example lines from `in.file`
-------------------
````
#Particle Number
particle_data number_particles 2
````

Explanation
------------
- `particle_data number_particles` equals number of particles
- Note that in case of a coupling to LIGGGHTS this value has to be equal to the total amount of particles in the simulation domain


Read next
-----------
 - [13_Coupling.md](13_Coupling.md)
