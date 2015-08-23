Coupling Model LIGGGHTS 
======================
This model is designed for coupling to the soft-sphere particle simulator 'LIGGGHTS'. The user needs to use  appropriate fixes in the LIGGGHTS simulation, in order to receive and supply data that is pulled and pushed to ParScale. For more details, see the example case supplied with the LIGGGHTS package.

Syntax
-------------------
````
    {
      "properties":
      {
        "verbose" : false,
        "pullEvery":  1
      }
    }
````

Explanation
------------
   - `verbose` enables a more detailed output
   - `pullEvery` sets the time steps where all important data is pulled and pushed to liggghts. Currently `pullEvery` has to be set to 1
     
     
   
IMPORTANT
------------
 - There must NOT be a 'control run' statement in the ParScale input script in case the user wants to run a simulation that is coupled to LIGGGHTS (the liggghts fix will not be found by ParScale when attempting to execute the control run).
     
Read next
-----------
 - [14_ModelEquation.md](14_ModelEquation.md)
 
 Related
-----------
 - [13_CouplingModel.md](13_CouplingModel.md) (for general information on coupling models)
