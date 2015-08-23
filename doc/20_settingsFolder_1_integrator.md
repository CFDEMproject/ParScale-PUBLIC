The settings/ folder - Integrator file
=======================================
For every simulations setting for the integrator have to be made in `settings/integrator.json` with the following syntax:

Syntax
-------------------

````
   {
    "integrator":       
    {
        "type":             "[String]",
        "linearSolver":     "[String]"
    },
    "parameters":
	{
      "absTol":     [Number],
      "relTol":     [Number],
      "mxsteps":    [Number],
      "maxord":     [Number],
      "deltaTInit": [Number],
      "deltaTMin":  [Number],
      "deltaTMax":  [Number],
      "maxNonLinearIterations": [Number]
    }
}
````

Example lines from `settings/integrator.json` 
-------------------

````    
    {
    "integrator":       
    {
        "type":             "CVODE",
        "linearSolver":     "CVDiag"
    },
    "parameters":
	{
      "absTol":     1e-6,
      "relTol":     1e-5,
      "mxsteps":    -1,
      "maxord":     3,
      "deltaTInit": 1e-6,
      "deltaTMin":  1e-12,
      "deltaTMax":  1.0,
      "maxNonLinearIterations": 2
    }
}
````

Note that at the moment only integrator of `type` `CVODE` and `linearSolver` of type `CVDiag` are available.

Read next
-----------
 - [30_0folder.md](30_0folder.md)