Installation Instructions
==================

Linux Environment
-----------------------

### Boost
In order to use the chemkin-II reader be sure you have the Boost library installed, e.g., on Ubuntu systems install 'boost-devel', 'libboost' and 'libboost-regex' packages (tested with version 1.54) with any common package manager. For a manual installation,

> visit [Libboost](http://www.boost.org/users/download/#live), download a verison of Libboost and follow the instruction

### HDF5
Be sure you have the HDF5 library installed, e.g., on OpenSUSE systems install all necessary hdf5 and hdf5-openmpi packages. For a manual installation, 

> visit [HDF5](http://www.hdfgroup.org/ftp/HDF5/current/src/unpacked/release_docs/INSTALL) and download a verison (*.gz) of HDF5

> `gunzip < hdf5-X.Y.Z.tar.gz | tar xf -`

> `cd hdf5-X.Y.Z`

> e.g. `./configure --prefix=/usr/local/hdf5 --enable-cxx`

> `make`

> `make install`

In case you do not wish to use HDF5 output, you can compile ParScale without HDF5 support. Simply use the `fedora_fpic_noHDF5` keyword (instead of `fedora_fpic`), as well as the `refresh*NoHDF5` scripts when compiling ParScale.

Also, it is useful to have a HDF5 viewer installed. Visit [this homepage](http://www.hdfgroup.org/products/java/release/download.html) to install 'hdfview'. Just download the installation script, and follow the installation instructions.

### Qt
Be sure you have the Qt library installed:

> Visit [QT Download Page](http://www.sysads.co.uk/2014/05/install-qt-5-3-ubuntu-14-04/) 

> Type (e.g. for 64 bit opperating system) `wget http://download.qt-project.org/official_releases/qt/5.3/5.3.0/qt-opensource-linux-x64-5.3.0.run`

> Type `chmod +x qt-opensource-linux-x64-5.3.0.run`

> Type `./qt-opensource-linux-x64-5.3.0.run`

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

Be sure you have correctly set the variables PASCAL_SRC_DIR, e.g., in your .bashrc you should have

>export PASCAL_LIGGGHTS_SRC_DIR=$HOME/LIGGGHTS/LIGGGHTS-TUG/src/

>export PASCAL_HDF5_DIR=/usr/local/hdf5

>export PASCAL_SRC_DIR=$HOME/LIGGGHTS/ParScale/src

>export PASCAL_EXAMPLE_DIR=$PASCAL_SRC_DIR/../examples

>export PASCAL_INST_DIR=$HOME/LIGGGHTS/ParScale/platforms/linux64

>export PASCAL_THIRDPARTY_DIR=$HOME/LIGGGHTS/ParScale/thirdParty/

>export PASCAL_SUNDIALS_DIR=$PASCAL_THIRDPARTY_DIR/sundials-2.5.0

>export PASCAL_QT5_DIR=$HOME/utilities/Qt/5.3/gcc_64

>export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PASCAL_QT5_DIR/lib:$PASCAL_THIRDPARTY_DIR/chemkinReader/src/:$PASCAL_HDF5_DIR/lib


General Hints
------------------------------------
- to build in parallel use "make -j numberOfProcessors" instead of just "make"

Building Third Party Libraries
------------------------------------
You need to build third party libraries before you build ParScale!
ParScale uses the following third party libraries:

- the [SUNDIALS package](http://computation.llnl.gov/casc/sundials). 
- a [ChemkinReader](https://github.com/lrm29/chemkinReader)

These libraries are included in the '/thirdParty/' directory of the ParScale source code package. To build these libraries, do a

>make thirdParty

To add more libraries of sub-packages to the compilation process, update the Makefile in this directory.

Building the Main Application
------------------------------------
ParScale will be build as an executable with just one instance of the PASCAL library. 
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
