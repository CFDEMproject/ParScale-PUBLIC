The settings/ folder - ParScale file
=======================================
For every simulations setting for the output produced by `ParScale` have to be made in `settings/parscale.json` with the following syntax:

Syntax
-------------------

````
{
    "input":
    {
        "writeJSON" : [Bool],
        "writeHDF5" : [Bool]
    }

}
````

Example lines from `settings/parscale.json` 
-------------------

````    
{
    "input":
    {
        "writeJSON" : true,
        "writeHDF5" : false
    }

}
````

The output will be printed in the specified format into the time folders. Not that most of the post processing and plotting scripts for `ParScale` require output in `.json` format.

Read next
-----------
 - [30_0folder.md](30_0folder.md)