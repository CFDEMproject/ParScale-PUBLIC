#!/bin/bash

#===================================================================#
# sytsem settings test routine for ParScale project 
# Christoph Goniva, March 2015
#===================================================================#

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#- system settings
clear
echo "*************************"
echo "ParScale system settings:"
echo "*************************"


echo ""
echo "checking if paths are set correctly..."
checkDirComment "$PASCAL_LIGGGHTS_SRC_DIR" '$PASCAL_LIGGGHTS_SRC_DIR' "yes"
checkDirComment "$PASCAL_HDF5_DIR" '$PASCAL_HDF5_DIR' "no"
checkDirComment "$PASCAL_SRC_DIR" '$PASCAL_SRC_DIR' "yes"
checkDirComment "$PASCAL_EXAMPLE_DIR" '$PASCAL_EXAMPLE_DIR' "no"
checkDirComment "$PASCAL_INST_DIR" '$PASCAL_INST_DIR' "yes"
checkDirComment "$PASCAL_THIRDPARTY_DIR" '$PASCAL_THIRDPARTY_DIR' "yes"
checkDirComment "$PASCAL_SUNDIALS_DIR" '$PASCAL_SUNDIALS_DIR' "yes"
checkDirComment "$PASCAL_QT5_DIR" '$PASCAL_QT5_DIR' "yes"
checkDirComment "$CFDEM_POEMSLIB_PATH" '$CFDEM_POEMSLIB_PATH' "yes"
echo ""

echo "echoing extra includes..."
echo PASCAL_EXTRA_INCLUDES: $PASCAL_EXTRA_INCLUDES

echo "echoing extra dirs amd libs..."
echo PASCAL_EXTRA_LIBS: $PASCAL_EXTRA_LIBS

echo ""
echo "Will use PASCAL_LIB_NAME=$PASCAL_LIB_NAME"

echo ""
exec="$PASCAL_SRC_DIR/$PASCAL_LIB_NAME"
if [ -f "$exec" ];then
    echo "The executable PASCAL_SRC_DIR/PASCAL_LIB_NAME = $PASCAL_SRC_DIR/$PASCAL_LIB_NAME exists."
else
    echo "No executable PASCAL_SRC_DIR/PASCAL_LIB_NAME found, you need to compile first."
fi

echo ""
echo "checking boost installation (only works for redhat, and if application dpkg is installed):"
if [ -f /etc/redhat-release ]; then
    echo "boost:"
    rpm -qa | grep boost
fi

if [ -f /etc/lsb-release ]; then
    echo "libboost-dev:"
    dpkg -s libboost-dev | grep 'Version'
    echo "libboost-regex-dev:"
    dpkg -s libboost-regex-dev | grep 'Version'
fi
