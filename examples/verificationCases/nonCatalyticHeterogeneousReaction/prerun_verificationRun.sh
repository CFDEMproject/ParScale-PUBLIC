#!/bin/bash

#Clean case
./Allclean

#Run the case
cp 0/speciesA_zero.json 0/speciesA.json
cp settings/chemistry_single_reaction_withSolidConsumption.json settings/chemistry_single_reaction.json
