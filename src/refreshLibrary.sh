#!/bin/sh

#***************************************
# Compiles ParScale LIBRARY
#***************************************

#copy coupling files from package to LIGGGHTS/src
cd $PASCAL_LIGGGHTS_SRC_DIR/PASCAL
chmod 777 Install.sh
./Install.sh 0
./Install.sh 1

#clean and remove all copulings to 
cd $PASCAL_SRC_DIR/COUPLING_LIGGGHTS
chmod 777 Install.sh
./Install.sh 0
./Install.sh 1
cd $PASCAL_SRC_DIR

make clean-all


#Compile
make thirdParty -j 4
make -f makelib
make -f Makefile.lib fedora_fpic -j 4
