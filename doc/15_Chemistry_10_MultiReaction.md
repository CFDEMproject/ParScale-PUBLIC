Chemistry MultiReaction
======================
The chemistry option `MultiReaction` enables the chemistry multi reacting model. Since ParScale is build for multi reaction problems, input in CHEMKIN-II format is required. 

Syntax
-------------------

````
    #Chemistry
    modelchemistry <modelName> reaction
````

Example lines from `in.file`
------------------------------

````
    #Chemistry
    modelchemistry MultiReaction reaction
```` 

 Explanation
----------------
- Files for reaction, thermophysics and transport are required which have to be named `chemkin.inp`, `therm.dat` and `tran.dat` and have to be placed in the running folder where the `in.file` is located
- Note that only reaction data is currently used for calculations
- Detailed output of the proceeded files are placed in the ParScale folder (`reactionsParsed`, `speciesParsed`)
- A detailed documentation of underlying models can be found in `doc/pdf/2_modelChem/`
- **For every species in `chemkin.inp` under SPECIES a modelEqn of type `species` has to be specified named after the reacting species (e.g. speciesOH for species of OH)**
- A example case including `chemkin.inp`, `therm.dat` and `tran.dat` can be found in `/home/tforg/ParScale/ParScale/examples/testCases/chemistryReader`


Read next
-----------
 - [16_EquationProperties.md](16_EquationProperties.md)
