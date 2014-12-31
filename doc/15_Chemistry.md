Chemistry Command
======================
The chemistry option enables the chemistry model. Since ParScale is build for multi reaction problems, input in CHEMKIN-II format is required. But also a simpler chemistry is available (`SingleReaction`). 
Note that if no chemistry model (e.g. non-reacting particles) is wanted the user just leaves out the `modelchemistry` command.

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
- `MultiReaction` and `SingleReaction` are the current available chemistry models (see separate *.md files)
- each modelchemistry requires the specification of a so-called "grain model". Grain models need to be specified in the file settings/chemistry_grainmodel.json. In case you do not wish to use a grain model, simply specify `"type" : "none"`. Also, you need to specify an entry for `chemistry_grainmodel` in the `settings/verbose.json` file in order to control the printing of debug messages to the screen.
- for a list of available grain models, see the `inline double` function definitions in the `src/chemistry_grainmodel.h` file


Read next
-----------
 - [15_Chemistry_1_SingleReaction.md](15_Chemistry_1_SingleReaction.md)

