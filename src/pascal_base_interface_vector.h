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

#ifndef PASC_BASE_INTERFACE_VECTOR_H
#define PASC_BASE_INTERFACE_VECTOR_H

#include "pascal_base_interface.h"
#include "pascal.h"
#include "psctype.h"
#include "mpi.h"
#include "error_base.h"
#include <string>
#include <vector>
#include <map>

//#define VERBOSE //developer to activate this if needed

using namespace std;

namespace PASCAL_NS
{

class ParScaleBaseInterfaceVector
{
  friend class ParScale;

  public:

    ParScaleBaseInterfaceVector()
    {
        error_ = 0;
    }

    void init() const
    {
        // HOW PASCAL WORKS

        // -1 make settings (any commands in input script)
        // -2 init (this function, called immediately before a run is executed)
        // -3 run (via run command), including phyiscs, coupling, output
        // - can repeat 1-3 as often as needed

        // order of execution defined in pascal.cpp Constructor


#ifdef VERBOSE
        printf("**allocating ... \n");
#endif
        for (unsigned i=0; i < list_.size(); i++)
            list_[i]->allocate();

        // read all data; either pull from LIGGGHTS (precedence)
        // or from JSON files; performed on proc 0
        // coupling_.read() is executed before particleData_.read()
        // see constructor of ParScale class in pascal.cpp
#ifdef VERBOSE
        printf("**reading ... \n");
#endif
        for (unsigned i=0; i < list_.size(); i++)
            list_[i]->read();

        // restart data from previous run if applicable
        // performed on proc 0
#ifdef VERBOSE
        printf("**getting restart data ... \n");
#endif
        for (unsigned i=0; i < list_.size(); i++)
            list_[i]->restart();

        // bcast all properties/settings so all MPI procs have it
#ifdef VERBOSE
        printf("**bcasting ... \n");
#endif
        for (unsigned i=0; i < list_.size(); i++)
            list_[i]->bcast();

        // parallellize - each proc just keeps data he needs
#ifdef VERBOSE
        printf("**parallelizing ... \n");
#endif
        for (unsigned i=0; i < list_.size(); i++)
            list_[i]->parallelize();

        // do model/class specific init after everything is bcasted
#ifdef VERBOSE
        printf("**initializing ... \n");
#endif
        for (unsigned i=0; i < list_.size(); i++)
            list_[i]->init();

#ifdef VERBOSE
        printf("** done! \n");
#endif
    }

    void parse_command(string map_string,int narg,char const* const* arg) const
    {
//        printf("narg %d, arg[0] %s\n",narg,arg[0]);
        if(narg < 1)
            error_->throw_error_all(FLERR,"Not enough arguments for :",map_string.c_str());

        if (map_.find(map_string) != map_.end())
        {
            map_[map_string]->parse_command(narg,arg);
        }
        else
            error_->throw_error_all(FLERR,"INPUT SCRIPT PARSING: class name not found:",map_string.c_str());
    }

  private:

    template<typename U>
    U* add(U *ptr,string map_string)
    {
      list_.push_back(static_cast<ParScaleBaseInterface*>(ptr));
      if(dynamic_cast<ErrorBase*>(ptr))
        error_ = dynamic_cast<ErrorBase*>(ptr);
      map_[map_string] = list_.back();
      if(ptr != list_.back())
        error_->throw_error_all(FLERR,"BAD INTERNAL ERROR: ASSERTION FAILED:",map_string.c_str());
      return static_cast<U*>(list_.back());
      //return ptr;
    }

    vector<ParScaleBaseInterface*> list_;
    mutable map<string,ParScaleBaseInterface*> map_;
    ErrorBase *error_;
};

} //end PASCAL_NS

#endif
