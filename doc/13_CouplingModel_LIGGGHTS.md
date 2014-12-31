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
     
     
   
     
Read next
-----------
 - [14_Model_Equations.md](14_Model_Equations.md)
 
 Related
-----------
 - [13_CouplingModel.md](13_CouplingModel.md) (for general information on coupling models)
