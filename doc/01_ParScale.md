ParScale (Particle Scale)
======================

ParScale is a open source code to resolve transient intra particle properties such as temperature and species profiles under various conditions. It gives the user the possibility to couple the program to other open source codes, such as the DEM code LIGGGHTS or apply transient boundary conditions by default. For more information of the possibilities to change relevant input please refer to the following documentation.

Every ParScale case has well defined required input from the user and a basic structure in the input which. One running case need to have the following components:

* Input file (later referred as `in.file`)
* `0/` folder for initial conditions
* `settings/` folder for model properties for all models specified by the user

If the user has set up the case and provide all the required input ParScale is executed with

````
    ParScale < in.file
````  

where `ParScale` is the application name.

Read next
-----------
 - [10_in.file.md](10_in.file.md)






