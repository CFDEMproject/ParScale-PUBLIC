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

	Parts of the code were developped in the frame of the NanoSim project funded
	by the European Commission through FP7 Grant agreement no. 604656.
\*-----------------------------------------------------------------------------------*/
#ifdef COUPLING_MODEL_CLASS

CouplingModelStyle(liggghts, CouplingModelLiggghts)

#else


#ifndef PASC_COUPLING_MODEL_LIGGGHTS_H
#define PASC_COUPLING_MODEL_LIGGGHTS_H

#include "coupling_model.h"
#include "lammps.h"
#include "modify.h"
#include "fix_pascal_couple.h"

namespace PASCAL_NS
{

class CouplingModelLiggghts : public CouplingModel
{
    public:

      CouplingModelLiggghts(ParScale *ptr, const char *_name);

      void init();

      void pull_n_bodies(int &_nbody, int &_nbody_all);
      void pull_box(double *_boxlo,double *_boxhi,double *_sublo,double *_subhi);
      void pull_proc_info(int *_procgrid,int *_myloc,int (&_procneigh)[3][2]);
      int* get_external_map(int &length);


    private:

      // pull
      bool fill_container_from_coupling(class ContainerBase &container) const;

      // push
      bool dump_container_to_coupling(class ContainerBase &container) const;

      bool isInitialized_ ;

      LAMMPS_NS::LAMMPS*                    lmp_;
      mutable LAMMPS_NS::FixParScaleCouple*   fix_coupling_;
};

} //end PASCAL_NS

#endif

#endif
