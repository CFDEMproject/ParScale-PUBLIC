#!/bin/bash

#Do the plotting
currDir=$PWD
cd octave 
octave < plotMe_Stage1.m
octave < plotMe_Stage2.m
octave < plotMe_conversionVsTime.m

#Show plot
#okular *.png &

#Return to current directory
cd $currDir
