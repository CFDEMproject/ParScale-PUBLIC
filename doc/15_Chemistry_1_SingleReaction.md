Chemistry Single Reaction
======================
The chemistry single reaction enables the user to resolve a particle with just one chemical reaction. It also give the possibility to study that reaction under the assumption of isothermal particle.
Note that no model heat equation is needed to perform that calculation since the temperature is set constant. In case the user  defines a heat equation PartScale will abort with the corresponding error message.

Syntax 
-------------------

````
    {
        "reactants":
        {
             "names": ["[String]",["[String]"]
        },
    
    
        "reaction":
    	{
    	  "cMinimum"       :  [number],
          "isIsoThermal"   :  [true/false],
          "temperature"    :  [number],  
          "Arrhenius_A"    :  [number],
          "Arrhenius_beta" :  [number] ,
          "Arrhenius_E_A"  :  [number]
        },
    
        "[String]":
        {
             "reactionOrder": [number],
             "stoichiometry": [number]
        },
        "[String]":
        {
             "reactionOrder": [number],
             "stoichiometry": [number]
        }
    }
````

Example from `chemistry_single_reaction.json`
------------------------------

````
    {
        "reactants":
        {
             "names": ["A","B"]
        },
    
    
        "reaction":
    	{
    	  "cMinimum"       :  1e-3,
          "isIsoThermal"   :  false,
          "temperature"    :  400,  
          "Arrhenius_A"    :  1e0,
          "Arrhenius_beta" :  0 ,
          "Arrhenius_E_A"  :  2000
        },
    
        "A":
        {
             "reactionOrder": 1,
             "stoichiometry": -1
        },
        "B":
        {
             "reactionOrder": 0,
             "stoichiometry": -1
        }
    }
```` 

 Explanation
----------------
- `reactants` are the names of the species taking part in the reaction
- `cMinimum` only relevant for 0-order reaction: this is the concentration below which a linear approximation will be used to stop the reaction (i.e., for very small concentrations a zero-order reaction must stop)
- `isIsoThermal` set the simulation to an isotermal state or not (`true` or `false` as input parameters)
- `temperature` is read in if `isIsoThermal` is set to `true`
- `Arrhenius_A`,`Arrhenius_beta`, `Arrhenius_E_A` sets the Arrhenuis constants for the single reaction
- **For every species defined under `names` a `reactionOrder` and `stoichiometry` has to be set**


Read next
-----------
 - [16_EquationProperties.md](16_EquationProperties.md)
