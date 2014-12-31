#!/bin/sh

#***************************************
# Compiles ParScale LIBRARY
#***************************************

#clean and remove all copulings to 
currDir=$PWD
cd COUPLING_LIGGGHTS
./Install.sh 0
./Install.sh 1
cd $currDir

make clean-all


#Compile
make thirdParty -j 4
make -f Makefile.lib fedora_fpic -j 4
cp *.h Obj_fedora_fpic
make -f Makefile.lib fedora_fpic -j 4
