Phase Change Model Equilibrium Sharp
==============================
Simple model for a liquid penetrating the particle where solid particles are dissolved. If liquid is reaching the saturation concentration solid particles are precipitated and distributed within the particle.

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
   "SaturationConcentration": 
    {
	  "Concentration": [number]
    }
}
````

Example from e.g. `liquidToSolidSolification.json`
------------------------------

````
{
	"involvedPhases": 
	{
		"phaseNames": 
		[
			"liquid",
			"solid"
		],
		"speciesNames":
		[
			"Suspension",
			"NanoParticles"
		]
	},
	"SaturationConcentration": 
    {
		"Concentration": 1.25
    }
}
```` 

 Explanation
----------------- 
- `phaseNames` specifies the "source" and "product" phase
- `speciesNames` refers the the user specified names of that phases
- `Concentration` within `SaturationConcentration` is setting the saturation concentration to a constant value


Read next
-----------
 -  [20_settingsFolder.md](20_settingsFolder.md)

