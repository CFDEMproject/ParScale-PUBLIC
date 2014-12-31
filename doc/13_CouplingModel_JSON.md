Coupling Model JSON
======================
This model is designed for offline usage, e.g., single-particle studies. Via this coupling model, the user can specify constant, or transient boundary conditions via a simulation run. Note, that the user can manipulate ANY data field (registered in the particle data class of ParScale).

Syntax
-------------------
````
    {
      "properties":
      {
        "verbose" : true,
        "pullEvery": [number]
      },
    
      "setBCs":
      {
          "[nameofsetting]" :
          {
            "dataField":  "[nameoffield]",
            "type"     :  "applyToAllParticles",
            "applyAt"  :  "[location]",
            "times"    :  [[number], [number], [number], [number]  ],
            "values"   :  [[number], [number], [number], [number]  ]
          },
          "[nameofsetting]" :
          {
            "dataField":  "[nameoffield]",
            "type"     :  "applyToAllParticles",
            "times"    :  [[number], [number]],
            "values"   :  [[number], [number]]
          }
      }
    }
````

Example lines from `coupling_json.json`
-----------------------------

````
    {
      "properties":
      {
        "verbose" : true,
        "pullEvery": 1
      },
    
      "setBCs":
      {
          "settingsHeat" :
          {
            "dataField":  "heat",
            "type"     :  "applyToAllParticles",
            "applyAt"  :  "fluid",
            "times"    :  [0,   0.005, 0.012, 0.1  ],
            "values"   :  [300, 400, 800, 900   ]
          },
          "settingsRadius" :
          {
            "dataField":  "radius",
            "type"     :  "applyToAllParticles",
            "times"    :  [0],
            "values"   :  [3e-3]
          }
      }
    }
````

Explanation
------------
- `verbose` enables a more detailed output
- `pullEvery` sets the time steps the wanted value is updated
- `[nameoffield]` is the setting name - more than one setting is possible with an arbitrary name (e.g. `myFancySetting`)
- `dataField` defines the particle data where the settings are applied (refer to the `<ModelEqnName>` the user sets in the Model Equation settings)
- `type` select the particles the settings are applied. Currently only `applyToAllParticles` is supported which will apply the settings to all particles
- `applyAt` defines at which intra-particle location the setting will be applied (only useful for array data, i.e., NOT required for scalar per-particle data). The user may specify 'surface' to make the setting directly at the surface. The setting 'fluid' will make the setting for the ambient fluid property
- `times` specifies the absolute time which in between the transient boundary conditions are calculated. `[0,10]` would indicate a calculation between 0 and 10 seconds
- `values` specifies the boundary values in between which are linearly interpolated
- E.g., in this example, at `time` = 0 s the boundary value of `heat` would be 300 K , at `time` = 0.005 s the boundary value (fluid temperature) of `heat` would be 400 K
     
     
   
     
Read next
-----------
 - [14_Model_Equations.md](14_Model_Equations.md)
 
 Related
-----------
 - [13_CouplingModel.md](13_CouplingModel.md) (for general information on coupling models)
