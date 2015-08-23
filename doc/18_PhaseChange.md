Phase Change Models
====================
Models in order to calculate phase change phenomena in `ParScale`.

Syntax
-------------------

````
    #Phase Change
    modelPhaseChange <modelName> <Name>
````

Example lines from `in.file`
------------------------------

````
    #Phase Change
    modelPhaseChange Evaporation liquidToGasEvaporation
```` 

 Explanation
----------------
 - `modelPhaseChange` indicates that there is a phase change model active.
 - `<modelName>` names the specific model. Current models are:
    * Evaporation
 - `<Name>` is the users name for the phase change after which the `.json` file in the `settings/` folder has to be named.

Read next
-----------
 -  [20_settingsFolder.md](20_settingsFolder.md)

Related
----------
 -  [18_PhaseChange_1_Evaporation.md](18_PhaseChange_1_Evaporation.md)