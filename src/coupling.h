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

#ifndef PASC_COUPLING_H
#define PASC_COUPLING_H

#include "error.h"
#include "stdio.h"
#include "coupling_base.h"
#include "coupling_model.h"
#include "pascal_base_accessible.h"
#include "pascal_base_interface.h"

namespace PASCAL_NS
{

class Coupling : public ParScaleBaseAccessible, public ParScaleBaseInterface, public CouplingBase
{
    public:

      Coupling(ParScale *ptr);
      ~Coupling();

      void read();

      void write();

      void bcast();

      void parallelize();

      void init();

      //TODO: CouplingModel
      // route functions below to CouplingModel

      virtual void parse_command(int narg,char const* const* arg);

      virtual void pull();
      virtual void pull_tags() {}; //TODO: Implement
      virtual void push(); 

      // check if active coupling exists
      bool external_code_in_control() const;

      // TODO extend in case of > 1 couplingmodel
      CouplingModel& couplingModel()
      { return *couplingModels_[0]; }

    private:

      mutable bool verbose_;

      // pull
      bool fill_container_from_coupling(class ContainerBase &container) const;

      // push
      bool dump_container_to_coupling(class ContainerBase &container) const;

      //Container for coupling models
      typedef CouplingModel *(*CouplingModelCreator)(ParScale*, char *name);
      std::map<std::string,CouplingModelCreator> *couplingModel_map_;

      template <typename T> static CouplingModel *couplingModel_creator(ParScale *ptr, char *name);

      vector<CouplingModel*> couplingModels_;

      bool isInitialized_;
};

} //end PASCAL_NS

#endif
