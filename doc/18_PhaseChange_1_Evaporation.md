Phase Change Model Evaporation
==============================
In order to correctly model evaporation in a multi-phase system a simple model for evaporation has been implemented which requires the following input parameters.

Syntax 
-------------------

````
    {
        "involvedPhases": {
                "phaseNames": [
                        "[String]",
                        "[String]"
                ],
                "speciesNames": [
                        "[String]",
                        "[String]"
                ]
        },
        "vaporPressureModel": {
                "evaporationRateConstant": [Number],
                "jacobiForDensePhase": 0.1,
                "pressureAmbient": [Number],
                "name": "[String]",
                "parameters": [[Number],[Number],[Number]],
                "source": "[String]",
                "comment": "[String]"
        },
        "enthalpyModel": {
                "name": "[String]",
                "value": [Number],
                "comment": "[String]"
        }
}
````

Example from e.g. `liquidToGasEvaporation.json`
------------------------------

````
    {
        "involvedPhases": {
                "phaseNames": [
                        "liquid",
                        "gas"
                ],
                "speciesNames": [
                        "liquidSpecies",
                        "species"
                ]
        },
        "vaporPressureModel": {
                "evaporationRateConstant": 1.0,
                "jacobiForDensePhase": 0.1
                "pressureAmbient": 100000,
                "name": "antoine",
                "parameters": [8.20417,1642.89,-42.85],
                "source": "",
                "comment": "gives 0.172 10e-5 [Pa] at 318.15 - works!model must return vapor pressure in Pa, temperature is in Kelvin"
        },
        "enthalpyModel": {
                "name": "constant",
                "value": 40500000,
                "comment": "Ethanol (65°C),model must return heat of evaporation in J/kmol and at 273 K"
        }
}
```` 

 Explanation
----------------- 
- `phaseNames` specifies the "source" and "product" phase
- `speciesNames` refers the the user specified names of that phases
- `evaporationRateConstant`, `jacobiForDensePhase`, `pressureAmbient` are parameters which refer to the physical model (see `pdf/1_modelEqn/00_base_model_eqn.pdf` for equations)
- `name` names the vapor pressure model used for the calculation of the pressure. The following models are available:
    * `antoine`: Antoine pressure model
- `parameters` are the parameters used for the approximation. Note that [DDBST](http://ddbonline.ddbst.de/AntoineCalculation/AntoineCalculationCGI.exe?) offers model parameters for a variety of liquids. Be careful to choose the last parameter (`C`) such as the equation is valid for a temperature in [K]
- `source` and `comment` gives the user space to place sources and add comments for other users or later usage 
- `name` of the `enthalpyModel` relates to the available model for the liquids enthalpy. The following models are implemented:
    * `constant`: constant enthalpy is assumed
- `comment` provides space for comments


Read next
-----------
 -  [20_settingsFolder.md](20_settingsFolder.md)

