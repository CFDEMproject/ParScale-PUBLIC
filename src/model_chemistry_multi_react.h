/*------------------------------------------------------------------------------------*\

                                      /$$$$$$                      /$$          
                                     /$$__  $$                    | $$          
        /$$$$$$   /$$$$$$   /$$$$$$ | $$  \__/  /$$$$$$$  /$$$$$$ | $$  /$$$$$$ 
       /$$__  $$ |____  $$ /$$__  $$|  $$$$$$  /$$_____/ |____  $$| $$ /$$__  $$
      | $$  \ $$  /$$$$$$$| $$  \__/ \____  $$| $$        /$$$$$$$| $$| $$$$$$$$
      | $$  | $$ /$$__  $$| $$       /$$  \ $$| $$       /$$__  $$| $$| $$_____/
      | $$$$$$$/|  $$$$$$$| $$      |  $$$$$$/|  $$$$$$$|  $$$$$$$| $$|  $$$$$$$
      | $$____/  \_______/|__/       \______/  \_______/ \_______/|__/ \_______/
      | $$                                                                      
      | $$                                                                      
      |__/        A Compilation of Particle Scale Models

   Copyright (C): 2014 DCS Computing GmbH (www.dcs-computing.com), Linz, Austria
                  2014 Graz University of Technology (ippt.tugraz.at), Graz, Austria
---------------------------------------------------------------------------------------
License
    ParScale is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with ParScale. If not, see <http://www.gnu.org/licenses/lgpl.html>.

	This code is designed to simulate transport processes (e.g., for heat and
	mass) within porous and no-porous particles, eventually undergoing
	chemical reactions.

	Parts of the code were developed in the frame of the NanoSim project funded
	by the European Commission through FP7 Grant agreement no. 604656.
\*-----------------------------------------------------------------------------------*/


//This is a class called if model_chemistry detects multi reaction (in CHEMKIN format) in settings/model_chemistry.json
//Holds all the calculations done with CHEMKIN format
//read in to be done by chemistry_reader_chemkin

#ifdef MODEL_CHEMISTRY_CLASS

ModelChemistryStyle(MultiReaction, ModelChemistryMultiReact)

#else

#ifndef PASC_CHEM_MULTI_REACT_H
#define PASC_CHEM_MULTI_REACT_H

#include "model_chemistry.h"
#include "model_base.h"
#include "model_container.h"
#include "error.h"
#include "pascal_base_accessible.h"
#include "chemistry_reaction_single.h"
#include <iostream>
#include "particle_data.h"
#include "model_eqn_container.h"
#include "chemistry_grainmodel.h"

namespace PASCAL_NS
{

class ModelChemistryMultiReact : public ModelChemistry
{
    public:

      ModelChemistryMultiReact(ParScale *ptr, char *name);
   
       
    private:
};

} //end PASCAL_NS

#endif

#endif

