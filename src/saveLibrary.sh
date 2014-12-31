#!/bin/sh

#***************************************
# Compiles ParScale LIBRARY
#***************************************

#clean and remove all copulings to 
currDir=$PWD
cp *liggghts* COUPLING_LIGGGHTS
cd COUPLING_LIGGGHTS
./Install.sh 0
cd $currDir

make clean-all

