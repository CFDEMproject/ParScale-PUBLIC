The settings/ folder
===============
For every property of defined in the `in.file`, a corresponding .json file in the `settings/` folder has to exist which specifies the value of the property.
So a model called `model_heatCapacity_solid` would result in a file called `model_heatCapacity_solid.json` with the following syntax:

Syntax
-------------------

````
    {
        "properties":
    	{
          "model" :  "const",
    	  "const_value" :  [number]
        }
    	
    }
````

Example lines from `settings/` (e.g. `model_heatCapacity_solid.json`)
-------------------

````    
    {
        "properties":
    	{
          "model" :  "const",
    	  "const_value" :  390.2
        }
    	
    }
````

 Explanation
----------------
 - currently only constant models for thermophysical properties are implemented so everything except the value of the property has to remain as shown
 - in case you are missing a required property file have a look through the examples and error messages displayed while executing your simulation 

Read next
-----------
 - [30_0Folder.md](30_0Folder.md)

Related
----------
 - [15_Chemistry.md](15_Chemistry.md)
 - [14_ModelEquation.md](14_ModelEquation.md)
 - [16_EquationProperties.md](16_EquationProperties.md)
 - [20_settingsFolder_1_integrator.md](20_settingsFolder_1_integrator.md)
 - [20_settingsFolder_2_parscale.md](20_settingsFolder_2_parscale.md)

