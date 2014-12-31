Model Equation Command
======================
The Model Equation option let the user define the equations in order to resolve the related property inside the particle. 

Syntax
-------------------

````
    #Model Equations
    modelEqn <ModelEqnType>  <ModelEqnName>     BC0 [number]  BC1 [number]
````

Example lines from `in.file`
-----------------------------

````
   #Model Equations
   modelEqn 1DSpherical     heat     BC0 1   BC1 2
   modelEqn ShrinkingCore   species          BC1 2

```` 
 Explanation
----------------
 - `modelEqn` indicates that a modelEqn option is following
 - For the currently implemented equation types, i.e., '<ModelEqnType>' see related documentation.
 - '<ModelEqnName>' allows the user to specify the name of an equation, which represents the name of the variable the model equation is solving for.
   
 
Read next
-----------
 - ***15_Chemistry.md***

Related
----------
 - [14_ModelEquation_1_1DSpherical.md](14_ModelEquation_1_1DSpherical.md)
 - [14_ModelEquation_2_ShrinkingCore.md](14_ModelEquation_2_ShrinkingCore.md)
