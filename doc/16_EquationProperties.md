Equation Properties
======================

Depending on the `<ModelEqnType>` model properties have to be defined and named with respect to the `<ModelEqnName>`. IN addition to that 


Syntax
-------------------

````
    #Heat properties
    model propertiesThermo <propertyName>
````

Example lines from `in.file`
-----------------------------

````
    #Heat properties
    model propertiesThermo heatThermalConductivity_solid 
    model propertiesThermo heatCapacity_solid 
    model propertiesThermo heatDensity_solid 
    model propertiesThermo heatTransferCoeff
    model propertiesThermo heatglobal_properties
```` 

 Explanation
----------------
- a `1DSpherical` with name `heat` would require the above input
- in case you are missing a required property have a look through the examples and error messages displayed while executing your simulation 
- a `1DSpherical` with name `speciesH2` would require the example input shown below

```` 
    model propertiesThermo speciesH2Diffusivity 
    model propertiesThermo speciesH2TransferCoeff 
    model propertiesThermo speciesH2global_properties
```` 

Read next
-----------
 - [17_Control.md](17_Control.md)

Related
----------
 - [15_Chemistry.md](15_Chemistry.md)
 - [14_ModelEquation.md](14_ModelEquation.md)



