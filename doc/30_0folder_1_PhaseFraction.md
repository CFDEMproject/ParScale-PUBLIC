The 0/ folder - Phase Fractions
===============================
Since `ParScale` is mainly designed for porous particles with three phases, phase fractions for every occurring phase have to be set. This requires the files `0/gasPhaseFraction.json` and/or `0/liquidPhaseFraction.json`. The syntax looks similar for e.g. `heat.json`. For particle at every grid point a phase fraction has to be specified.  

Syntax
-------------------

````
    {
        "name": "`gasPhaseFraction/liquidPhaseFraction`",
        "data" :
    {
        "1": [ [number], [number], [number], [number], [number], [number],[BCvalue]],
        "2": [ [number], [number], [number], [number], [number], [number],[BCvalue]]
    }
    }
 
````

Example lines from `0/gasPhaseFraction.json`
-----------------------------

````
    {
        "name": "gasPhaseFraction",
        "data" :
        {
            "1": [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5],
            "2": [0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3]
        }
    }
    
Read next
-----------
 - [40_testHarnessConfigFiles.md](40_testHarnessConfigFiles.md)
