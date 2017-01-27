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

     Copyright (C): 2012 - 2014 DCS Computing GmbH (www.dcs-computing.com), Linz, Austria
                    2012 - 2014 Department of Particulate Flow Modelling, JKU Linz
                              (www.jku.at/pfm), Linz, Austria

   This file was originally part of LIGGGHTS (www.cfdem.com), and is now re-distributed
   under LGPL as part of ParScale with the permission of the copyright holders
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

#include "custom_value_tracker.h"
#include "input.h"
#include "output.h"
#include "particle_data.h"
#include "coupling.h"
#include "comm.h"
#include "mpi_pascal.h"

using namespace PASCAL_NS;

  /* ----------------------------------------------------------------------
   constructor, destructor
  ------------------------------------------------------------------------- */

  CustomValueTracker::CustomValueTracker(ParScale *ptr, ParticleData *pdata)
   : ParScaleBase(ptr),
     nbody_(0),
     nbody_all_(0),
     owner_(*pdata),
     capacityElement_(0)
  {
  }

  CustomValueTracker::~CustomValueTracker()
  {
  }

  /* ----------------------------------------------------------------------
   memory management
  ------------------------------------------------------------------------- */

  int CustomValueTracker::getCapacity()
  {
    return capacityElement_;
  }

  /* ----------------------------------------------------------------------
   remove property
  ------------------------------------------------------------------------- */

  void CustomValueTracker::removeElementProperty(const char *_id)
  {
     elementProperties_.remove(_id);
  }

  void CustomValueTracker::removeGlobalProperty(const char *_id)
  {
     globalProperties_.remove(_id);
     globalProperties_orig_.remove(_id);
  }

  /* ----------------------------------------------------------------------
   initial setting of nbody
  ------------------------------------------------------------------------- */

  void CustomValueTracker::set_n_body(int _nbody,int _nbody_all)
  {

    nbody_     = _nbody;
    nbody_all_ = _nbody_all;

#if 0
    printf("Will track %d particles (%d global).\n \n",
           nbody_,nbody_all_);
#endif
  }

  /* ----------------------------------------------------------------------
   initial allocation
  ------------------------------------------------------------------------- */

  void CustomValueTracker::allocate()
  {

    for(int iP = 0; iP < nbody_; iP++)
        addZeroElement();

    char msgstr[500];
    sprintf(msgstr,
            "CustomValueTracker::allocate added %d particles (of %d global) to the simulation...\n",
            nbody_,nbody_all_);
    output().write_screen_all(msgstr);

  }

  /* ----------------------------------------------------------------------
   initial allocation
  ------------------------------------------------------------------------- */

  void CustomValueTracker::recalc_nbody_all(bool errflag)
  {
    int nbody_all_old = nbody_all_;
    MPI_Sum_Scalar(nbody_,nbody_all_,comm().world());
    if(errflag && nbody_all_ != nbody_all_old)
        error().throw_error_all(FLERR,"Particles were lost");
  }

  /* ----------------------------------------------------------------------
   store current values of global properties as orig
  ------------------------------------------------------------------------- */

  void CustomValueTracker::storeOrig()
  {
      //NP this handles owned and ghost elements
      globalProperties_.storeOrig(globalProperties_orig_);
      /*NL*/ //fprintf(screen,"storeOrig() called \n");
      /*NL*///  error->all(FLERR,"Internal error");
  }

  /* ----------------------------------------------------------------------
   reset global properties to orig
  ------------------------------------------------------------------------- */

  void CustomValueTracker::resetToOrig()
  {
      //NP this handles owned and ghost elements
      globalProperties_.reset(globalProperties_orig_);
      /*NL*/ //fprintf(screen,"resetToOrig() called \n");
      /*NL*/ //error->all(FLERR,"Internal resetToOrig called");
  }

  /* ----------------------------------------------------------------------
   rotate all properties, applies to vector and multivector only
  ------------------------------------------------------------------------- */

  void CustomValueTracker::rotate(double *totalQ,double *dQ)
  {
      /*NL*/ //printVec4D(screen,"totalQ",totalQ);
      /*NL*/ //printVec4D(screen,"dQ",dQ);

      //NP this handles owned and ghost elements
      elementProperties_.rotate(dQ);
      globalProperties_.rotate(totalQ);
  }

  void CustomValueTracker::rotate(double *dQ)
  {
      //NP this handles owned and ghost elements
      elementProperties_.rotate(dQ);
      globalProperties_.rotate(dQ);
  }

  /* ----------------------------------------------------------------------
   scale all properties, applies to vectors and multivectors only
  ------------------------------------------------------------------------- */

  void CustomValueTracker::scale(double factor)
  {
      //NP this handles owned and ghost elements
      elementProperties_.scale(factor);
      globalProperties_.scale(factor);
  }

  /* ----------------------------------------------------------------------
   move all properties
  ------------------------------------------------------------------------- */

  void CustomValueTracker::move(double *vecTotal, double *vecIncremental)
  {
      //NP this handles owned and ghost elements
      elementProperties_.move(vecIncremental);
      globalProperties_.move(vecTotal);
  }

  void CustomValueTracker::move(double *vecIncremental)
  {
      //NP this handles owned and ghost elements
      elementProperties_.move(vecIncremental);
      globalProperties_.move(vecIncremental);
  }

  /* ----------------------------------------------------------------------
  re-sort data stored in element properties
  ------------------------------------------------------------------------- */

  void CustomValueTracker::sortPropsByExtMap(int *_id, int _nlocal,
                                             int &_len_id,
                                             int *_map, int _len_map,
                                             bool verbose, int me)
  {
    elementProperties_.sortPropsByExtMap( _id, _nlocal,
                                          _len_id,
                                          _map, _len_map,&error(),
                                          verbose, me
                                        );
  }

  /* ----------------------------------------------------------------------
   clear reverse properties, i.e. reset all of them to 0
  ------------------------------------------------------------------------- */

  void CustomValueTracker::clearReverse(OperationProperties &op)
  {
      //NP this handles owned and ghost elements
      elementProperties_.clearReverse(op);
  }

  /* ----------------------------------------------------------------------
   read from JSON files
  ------------------------------------------------------------------------- */

  void CustomValueTracker::read(OperationProperties &op)
  {
      elementProperties_.read(op,static_cast<InputBase const*>(&input()));
      globalProperties_.read(op,static_cast<InputBase const*>(&input()));
  }

  /* ----------------------------------------------------------------------
   pull from LIGGGHTS
  ------------------------------------------------------------------------- */

  void CustomValueTracker::pull(OperationProperties &op)
  {
      elementProperties_.pull(op,static_cast<CouplingBase const*>(&coupling()));
      globalProperties_.pull(op,static_cast<CouplingBase const*>(&coupling()));
  }

  /* ----------------------------------------------------------------------
   write JSON files
  ------------------------------------------------------------------------- */

  void CustomValueTracker::write(OperationProperties &op)
  {
      elementProperties_.write(op,static_cast<InputBase const*>(&input()));
      globalProperties_.write(op,static_cast<InputBase const*>(&input()));
  }

  /* ----------------------------------------------------------------------
   push to LIGGGHTS
  ------------------------------------------------------------------------- */

  void CustomValueTracker::push(OperationProperties &op)
  {
      elementProperties_.push(op,static_cast<CouplingBase const*>(&coupling()));
      globalProperties_.push(op,static_cast<CouplingBase const*>(&coupling()));
  }
