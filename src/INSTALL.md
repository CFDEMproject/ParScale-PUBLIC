Installation Instructions
==================

Linux Environment
-----------------------

### Boost
In order to use the chemkin-II reader be sure you have the Boost library installed, e.g., on Ubuntu systems install 'boost-devel', 'libboost' and 'libboost-regex' packages (tested with version 1.54) with any common package manager. For a manual installation,

> visit [Libboost](http://www.boost.org/users/download/#live), download a verison of Libboost and follow the instruction

After installation, be sure you can link the shared object 'lboost_regex', e.g., by typing 'gcc -lboost_regex'. Also, make sure the libboost_regex.so file is availabe on each cluster node (ask you system admin, or make a copy of the library on your home directory). In some cases, the user may need to set the variable 'PASCAL_EXTRA_INCLUDES' and 'PASCAL_EXTRA_LIBS' to successfully link boost installation files during the compilation process (see instructions below).

### HDF5
Be sure you have the HDF5 library installed, e.g., on OpenSUSE systems install all necessary hdf5 and hdf5-openmpi packages. 

WARNING: Be sure you have HDF5 version 1.8.15 or higher installed!

For a manual installation, 

> visit [HDF5](http://www.hdfgroup.org/ftp/HDF5/current/src/unpacked/release_docs/INSTALL) and download a verison (*.gz) of HDF5

> `gunzip < hdf5-X.Y.Z.tar.gz | tar xf -`

> `cd hdf5-X.Y.Z`

> e.g. `./configure --prefix=/usr/local/hdf5 --enable-cxx`

> `make`

> `make install`

After installation of hdf5, make sure the directories $PASCAL_HDF5_DIR/include and  $PASCAL_HDF5_DIR/lib or  ($PASCAL_HDF5_DIR/lib64) exist and contain include files and all relevant shared libraries.

In case you DO NOT wish to use HDF5 output, you can compile ParScale without HDF5 support. Simply use the `fedora_fpic_noHDF5` keyword (instead of `fedora_fpic`), as well as the `refresh*NoHDF5` scripts when compiling ParScale. Also, set the `PASCAL_HDF5_DIR` to empty in your .bashrc, i.e., do

> `export PASCAL_HDF5_DIR=""`

in your .bashrc file to let your operating system know that ParScale does not need HDF5 libraries during the linking process.

Also, it is useful to have a HDF5 viewer installed. Visit [this homepage](http://www.hdfgroup.org/products/java/release/download.html) to install 'hdfview'. Just download the installation script, and follow the installation instructions.

### Qt
Be sure you have the Qt library installed:

> Visit [QT Download Page](http://download.qt.io/official_releases/qt/5.8/5.8.0/) 

> Type (e.g. for 64 bit opperating system) `wget http://download.qt.io/official_releases/qt/5.8/5.8.0/qt-opensource-linux-x64-5.8.0.run`

> Type `chmod +x qt-opensource-*.run`

> Type `./qt-opensource-*.run`

WARNING: if you use a QT version 5.6 >, you may require C++11 support. To activate this support when compiling, the user may set the environment variable "PASCAL_C11_STD" appropriately. For example, adding 

> export PASCAL_C11_STD=-std=c++11

to the .bashrc file will enable the C++11 support when compiling with gcc/mpic++.

### Octave & JSONLAB
We recommend using octave version 3.8.x or later for correct display and printing of result graphs.

Be sure you have correctly set up Octave (including `JSONLAB`) for post processing

> Visit [jsonlab Download Page](http://sourceforge.net/projects/iso2mesh/files/jsonlab/) 

> Download and unpack latest Version (here we assume that you save it to `'/home/username/utilities/jsonlab'`)

> Go to your $HOME folder and type `touch .octaverc`

> Edit `.octaverc` (e.g. `vim .octaverc`) and add the path to jsonlab folder by typing e.g. `addpath('/home/username/utilities/jsonlab')`

> Save `.octaverc`, close and countercheck the added path by typing `octave`

> In Octave type `path` to check all loaded paths - `/home/username/utilities/jsonlab` should appear on top


### Environmental Variables

Be sure you have correctly set the variables ParScale requires for compilation and running, i.e., in your .bashrc you should have (note, the setting for PASCAL_C11_STD may depend on your QT version, since some require C++11 coding standards)

>export PASCAL_LIGGGHTS_SRC_DIR=$HOME/LIGGGHTS/LIGGGHTS-TUG/src/

>export PASCAL_HDF5_DIR=/usr/local/hdf5

>export PASCAL_SRC_DIR=$HOME/LIGGGHTS/ParScale/src

>export PASCAL_QT5_DIR=$HOME/utilities/Qt/5.3/gcc_64

>export PASCAL_C11_STD=-std=c++11

>export PASCAL_LIB_NAME=pasc_fedora_fpic

>export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PASCAL_QT5_DIR/lib:$PASCAL_HDF5_DIR/lib64:$PASCAL_THIRDPARTY_DIR/chemkinReader/src/

The user may also want to use a bashrc file to define some environmental variables, see here:
>. $PASCAL_SRC_DIR/../etc/bashrc

The user can also set 'PASCAL_EXTRA_INCLUDES' and 'PASCAL_EXTRA_LIBS' in the .bashrc file to add include path and libraries needed to compile ParScale (for example, this is useful if the user is unable to install the boost libraries as root). Please also see the instructions in the file '$PASCAL_SRC_DIR/../etc/bashrc'. Note that the 'PASCAL_EXTRA_INCLUDES' and 'PASCAL_EXTRA_LIBS' variables are accessed AS IS in the make scripts, so the user MUST INCLUDE the correct symbols in the variable definition. For example, a valid setting in the .bashrc file is

> export PASCAL_EXTRA_INCLUDES="-I$PASCAL_QT5_DIR/include/qt5"
> export PASCAL_EXTRA_LIBS="-L$HOME/lib"

The user should run the $PASCAL_SRC_DIR/../etc/parScaleSystemTest.sh and check that all path are set correctly.

General Hints
------------------------------------
- to build in parallel use "make -j numberOfProcessors" instead of just "make"

Building Third Party Libraries
------------------------------------
Source code for the required third party libraries is included in the '/thirdParty/' directory of the ParScale source code package. In summary, ParScale uses the following third party libraries:

- the [SUNDIALS package](http://computation.llnl.gov/casc/sundials). 
- a [ChemkinReader](https://github.com/lrm29/chemkinReader)

You need to build these third party libraries before you build ParScale using ParScale's make system. In case you run one of the shell scripts in the src directory (e.g., 'refresh', see below), third party libraries will be already build, so you typically do not need to worry about building these libraries! To build these libraries manually, do a

>make thirdParty

To add more libraries of sub-packages to the compilation process, update the Makefile in this directory.

Building the Main Application
------------------------------------
ParScale will be built as an executable with just one instance of the PASCAL library. 
Currently, only one linux architecture is supported.  

The script 

>src/refresh.sh

will do a make clean-all first, and remove all couplings to external codes (e.g., LIGGGHTS). Then it will compile third-party software, as well as ParScale.

To compile manually, simply do:

>make fedora_fpic

Note, a

>make clean-all 

is recommended to clean out pre-compiled files. Note: in case dependencies are missing, this might be due to the fact that all couplings to external codes have not been removed. Check the src/refresh.sh script in such a case!

The above procedure will generate an executable (called `pasc_*`) in the src directory. After building, you should generate a symbolic link to the compiled executable (e.g., in the terminal do: `ln -s /home/<userName>/ParScale/src/pasc_fedora_fpic /home/<userName>/bin/ParScale`). 

Building the Library
------------------------------------
ParScale will be build as library. Currently, only one linux architecture is supported. 

The script 

>src/refreshLibrary.sh

will install all couplings (e.g., for LIGGGHTS), clean the code, and compile the code as a library.

To compile manually, simply do:

>make makelib

>make -f Makefile.lib fedora_fpic
 
Note: In case you have modified the coupling routines (e.g., for LIGGGHTS), use the script 

> src/saveLibrary.sh

to save the modified code to the respective COUPLING_* directory.
 
Note, a

>make clean-all 

is recommended to clean out ALL (including third party libraries) pre-compiled files. 

Note, that there might exist additional scripts in the `src` directory that can be used to build ParScale without specific sub-libraries (e.g., without the HDF5 output module).
