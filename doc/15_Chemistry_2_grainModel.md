Grain Model 
======================
ParScale has the option to set a grain model for the chemistry (for detailed equations see doc/ and source code). Whenever a chemistry model is declared `chemistry_grainmodel.json` has to be added in the `settings/` folder.

Syntax
-------------------

````
    {
        "modelParameters":
        {
             "type"         : "[model_name]",
             "solidID"      : [number],
             "cSolidInit"   : [number],
             "psi0"         : [number],
             "alpha"        : [number]       
        }
    }
````

Example lines from `chemistry_grainmodel.json`
------------------------------

````
    {
        "modelParameters":
        {
             "type"         : "grainClassical",
             "solidID"      : 1,
             "cSolidInit"   : 1.0,
             "psi0"         : 0.75,
             "alpha"        : 0.5       
        }
    }
```` 

 Explanation
----------------
- ´type´ is the model name. Available are: `grainClassical`, `volumetric`, `randomPore` and is no grain model is wanted `none`
- `solidID` is related to the ID of the solid component and has to be set to the corresponding model equation ID
- `cSolidInit` is the initial concentration of the solid
- `psi0`, `alpha` are model parameters which are explained in the documentation (`doc/pdf/`)


Read next
-----------
 - [16_EquationProperties.md](16_EquationProperties.md)
