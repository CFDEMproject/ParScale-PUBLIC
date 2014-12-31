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

#ifndef PASC_CHEM_GRAINMODEL_H
#define PASC_CHEM_GRAINMODEL_H


#include "stdio.h"
#include "error.h"
#include "input_base.h"
#include "pascal_base_accessible.h"
#include "container.h"
#include "pascal_base_interface.h"
#include "particle_mesh.h"
#include "model_eqn.h"
#include "model_eqn_container.h"
#include "particle_data.h"
#include "input.h"

namespace PASCAL_NS
{

class ChemistryGrainModel : public ParScaleBaseAccessible, ParScaleBaseInterface 
{
  public:

      ChemistryGrainModel(ParScale *ptr);
      ~ChemistryGrainModel();

    void  init();

    int   solidID() const {return solidID_;};
    
    bool  active() {return active_;};

    double f(double cSolid) const 
    {
        if(verbose_)
            printf("GrainModel for X = %g: %g. \n",
                    max(0.0,1.0-cSolid/cSolidInit_),
                    (this->*FX)( max(0.0,1.0-cSolid/cSolidInit_) )
                  );
        
        return (this->*FX)( max(0.0,1.0-cSolid/cSolidInit_) );
    };
    
    void printStatus() const;

  private:

    bool            verbose_; 
    QJsonObject     global_properties_;
    QJsonObject     grainmodel_properties_;

    //pointer to the specific grain model used
    double    (ChemistryGrainModel::*FX)(double conversion) const; 
    
    //Model Parameters
    bool            active_;       //general switch
    int             solidID_;      //id of the solid species 
    double          cSolidInit_;
    double          psi0_;
    double          alpha_; 

    
    // - - - - - - -  - - - - - - -  - - - - - - -  - - - - - - -  - - - - - - - 
    //Individual Grain Models
    // - - - - - - -  - - - - - - -  - - - - - - -  - - - - - - -  - - - - - - - 
    //No model, type: 'none'
    inline double  fxNone(double conversion) const 
    {     return 1.0; }
         
    //Szekely et al. (1976), type: 'volumetric'
    inline double  fxVolumetric(double conversion) const 
    {     return 1.0-conversion; }

    //van den Aarsen (1985), type: 'grainClassical'
    inline double  fxGrainClassical(double conversion) const 
    {     return pow(1.0-conversion,0.666666666666666666667); }


    //Bhatia and Perlmutter (1980), type: 'randomPore'
    inline double  fxRandomPore(double conversion) const 
    {     return    (1.0 - conversion)
                  *sqrt( 1.0 -   psi0_
                                *log(1.0-conversion)
                       ); 
    }

    //General power law model (Melchiori and Canu, 2014), type: 'powerLaw'
    inline double  fxPowerLaw(double conversion) const 
    {     return pow(1.0-conversion, alpha_); }

};

} //end PASCAL_NS

#endif
