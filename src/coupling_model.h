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

#ifndef PASC_COUPLING_MODEL_H
#define PASC_COUPLING_MODEL_H

#include "coupling_base.h"
#include "pascal_base.h"
#include "qjson_includes.h"

namespace PASCAL_NS
{

class CouplingModel : public ParScaleBase, CouplingBase, ParScaleBaseInterface
{
    public:

      CouplingModel(ParScale *ptr, const char *_name);
      ~CouplingModel();

      virtual void    init() {};

      virtual void    read();

      // pull
      virtual bool fill_container_from_coupling(class ContainerBase &container) const = 0;

      // push
      virtual bool dump_container_to_coupling(class ContainerBase &container) const = 0;

      //Access/Setting functions
      const char* name()      const {return name_;}
      const char* scope()     const {return scope_;}

      void  setPushPull()      const;

      std::vector<int>* exchangeEventsLocalId() const
      {
            return exchangeEventsLocalId_;
      }
      std::vector<int>* exchangeEventsReceivingProcess() const
      {
            return exchangeEventsReceivingProcess_;
      }

      bool verbose() const {return verbose_;};

      virtual bool external_code_in_control() const
      { return true; }

      virtual void pull_n_bodies(int &_nbody, int &_nghost, int &_nbody_all)
      { _nbody = _nghost = _nbody_all = 0; }

      virtual void pull_box(double *_boxlo,double *_boxhi,double *_sublo,double *_subhi) {}
      virtual void pull_proc_info(int *_procgrid,int *_myloc,int (&_procneigh)[3][2]) {}
      virtual void pull_timeStepping_info(double &deltaT, int &neighAgo, int &timeStepFromRun) {}
      virtual int* get_external_map(int &length) {return NULL;};

    private:

      char *name_;
      char *scope_; // scope/name for JSON file

      //holds general properties of the coupling
      QJsonObject properties_;

    protected:

      mutable bool verbose_;
      mutable int  pullEvery_;
      mutable int  pullNext_;
      mutable bool hasPulled_;

      mutable int  pushEvery_;
      mutable int  pushNext_;
      mutable bool hasPushed_;

      std::vector<int> * exchangeEventsLocalId_;
      std::vector<int> * exchangeEventsReceivingProcess_;
};

} //end PASCAL_NS

#endif
