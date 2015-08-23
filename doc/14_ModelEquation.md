Model Equation Command
======================
The Model Equation option let the user define the equations in order to resolve the related property inside the particle. 

Syntax
-------------------

````
    #Model Equations
    modelEqn <ModelEqnType>  <ModelEqnName>   <State>  BC0 [number]  BC1 [number] #extra_options
````

Example lines from `in.file`
-----------------------------

````
   #Model Equations
   modelEqn 1DSpherical     heat           BC0 1   BC1 2 <extra_option>
   modelEqn ShrinkingCore   species gas            BC1 2 
   modelEqn 1DSpherical     species liquid         BC1 2 <extra_option>

```` 
 Explanation
----------------
 - `modelEqn` indicates that a modelEqn option is following
 - For the currently implemented equation types, i.e., `<ModelEqnType>` see related documentation.
 - `<ModelEqnName>` allows the user to specify the name of an equation, which represents the name of the variable the model equation is solving for.
 - `<State>` defines the physical state of every equation of `<ModelEqnType>` `species`. `heat` does not require any `<State>`
 - `<extra_option>` can be set in order to activate more options for each model equation. Available options are:
    * `updatePhaseFraction`: if `<ModelEqnType>` is `species` and a phase change model is specified the user can activate `updatePhaseFraction` in order to update the phase fraction. The concentration of the model equation will not be changed
    * `verbose`: active for more detailed "on-screen" information of the calculations/integrations
    * `solveConvectiveFlux: activate to solve the convective part of the model equation - see documentation 
    * `writeDebugContainers`: useful information for debugging and finding input errors
    * `averagePhaseFraction`: writes the averaged values for phase fractions to a file
    * `inactive`: useful to deactivate e.g. a heat equation in order to achieve isothermal behavior since this equation is not solved 
- See documentation in `doc/pdf/1_modelEqn/00_base_model_eqn.pdf` for more infos and underlying equations.

Read next
-----------
 - [15_Chemistry.md](15_Chemistry.md)

Related
----------
 - [14_ModelEquation_1_1DSpherical.md](14_ModelEquation_1_1DSpherical.md)
 - [14_ModelEquation_2_ShrinkingCore.md](14_ModelEquation_2_ShrinkingCore.md)
