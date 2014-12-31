DEVELOPERS GUIDE
==================

ParScale has been implemented in a fashion similar to LAMMPS and LIGGGHTS. This has been done to make it easier for developers familiar with these widely used packages to modify, improve, or extend ParScale. For developers not familiar with LAMMPS and/or LIGGGHTS, we recommend to study the [LIGGGHTS-PUBLIC Wiki] Pages (https://github.com/CFDEMproject/LIGGGHTS-PUBLIC/wiki).

ParScale uses two essential libraries, namely Qt 5.x and HDF5, in order to streamline the implementation. Also, ParScale relies on third-party libraries that are delivered in the `thirdParty` directory of the ParScale source code package.

JSON Input/Output
-----------------------
We rely on the QJson* classes available from the Qt Project:

> Qt 5.x  [Qt Downloads](http://qt-project.org/downloads)


Here is an introduction to [QJson stuff](http://qt-project.org/doc/qt-5/json.html#the-json-classes). You can verify if you are usign correct JSON files using [JSONLint](http://jsonlint.com/)

HDF5 OUTPUT
-----------------------
Every time folder contains a sub-folder where the HDF5 Data is written. It can easily be read by typing e.g.
> h5dump heat.h5`

Also, note that there are GUIs that can be used to view or modify HDF5 files, e.g., [hdfview](http://www.hdfgroup.org/products/java/release/download.html)
For more information see [HDF5 Group](http://www.hdfgroup.org/) and for post processing with octave visit [Octave-HDF5](http://hdfeos.org/software/octave.php).
