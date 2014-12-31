Coupling Model
======================
The Coupling Model option gives the user the possibility to let ParScale interact with other open source codes (currently only LIGGGHTS is supported) or to apply transient boundary conditions to the particles (e.g. rising temperatures). Note, that only one coupling model can be specified per simulation!

Syntax
-------------------
````
    `coupling <couplingName> `
````

Example lines from `in.file`

````
    #Coupling
    coupling none
````

Explanation
------------
- `coupling` specifies the selected coupling model

Currently supported coupling models are:

        - none
        - json
        - liggghts

Through **none**, no coupling model is activated. **json** activates a coupling model to `.json` input files in order to activate transient boundary conditions. **liggghts** activates the coupling model to the open source code LIGGGHTS. 
Every coupling model needs a settings file in the `settings/` folder which is named after the <couplingName> (e.g. `coupling_liggghts` for liggghts coupling). E.g. for a **json** coupling the `coupling_json.json` looks like this:
   
     
Read next
-----------
 - [14_Model_Equations.md](14_Model_Equations.md)
 
 Related
-----------
 - [40_LIGGGHTSCoupling.md](40_LIGGGHTSCoupling.md) (for more information on coupling models)
