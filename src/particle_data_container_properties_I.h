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

#ifndef PSC_PARTICLE_DATA_CONTAINER_PROPERTIES_I_H
#define PSC_PARTICLE_DATA_CONTAINER_PROPERTIES_I_H

#include "string.h"

  /* ----------------------------------------------------------------------
   decide if property is pushed or pulled at all
  ------------------------------------------------------------------------- */

  inline bool ParticleDataContainerProperties::decidePackUnpackOperation(OperationProperties &op)
  {
      // return true for manual communication, such as for node_, node_orig_
      // etc in MultiNodeMeshParallel
      if(COMM_TYPE_MANUAL == communicationType_)
        return true;

      //NP check for restart
      if(OPERATION_RESTART == op.operation())
      {
          if(RESTART_TYPE_YES == restartType_)
            return true;
          return false;
      }

      if(OPERATION_READ == op.operation())
      {
        if(do_read_)
            return true;
        return false;
      }

      if(OPERATION_PULL == op.operation())
      {
        if( needsPull() )
            return true;
        return false;
      }
      
      if(OPERATION_PUSH == op.operation())
      {
        if( needsPush() )
            return true;
        return false;
      }

      //NP communication in exchange() and borders() steps is always performed
      if(OPERATION_COMM_BORDERS == op.operation() ||
         OPERATION_COMM_EXCHANGE == op.operation() )
        return true;

      //NP only containers that are read in from file have to be broadcasted
      if(OPERATION_COMM_BCAST == op.operation())
      {
        if(needsBcast_)
            return true;
        return false;
      }

      //NP no fw or reverse comm for COMM_TYPE_NONE, but do exchange or borders
      if(COMM_TYPE_NONE == communicationType_)
        return false;

      if(OPERATION_COMM_REVERSE == op.operation() &&
         COMM_TYPE_REVERSE == communicationType_)
        return true;

      if(OPERATION_COMM_FORWARD == op.operation() &&
         COMM_TYPE_FORWARD == communicationType_)
        return true;

      if(OPERATION_COMM_FORWARD == op.operation() &&
         COMM_TYPE_FORWARD_FROM_FRAME == communicationType_)
      {
         if(op.scale() && !isScaleInvariant())
           return true;
         if(op.translate() && !isTranslationInvariant())
           return true;
         if(op.rotate() && !isRotationInvariant())
           return true;

         return false;
      }

      if(OPERATION_UNDEFINED== op.operation())
      {
        printf("ERROR: ParticleDataContainerProperties::decidePackUnpackOperation() has to handle an object with OPERATION_UNDEFINED. Will return 'false' and proceed.");
        return false;
      }

      printf("ERROR: illegal call to ParticleDataContainerProperties::decidePackUnpackOperation(). Will return 'false' and proceed.");

      // default
      return false;
  }

  /* ----------------------------------------------------------------------
   decide if operation performs data communication
  ------------------------------------------------------------------------- */

  inline bool ParticleDataContainerProperties::decideCommOperation(OperationProperties &op)
  {
      //NP have to decide at unpack if data is initialized with 0
      //NP or pulled from buffer

      //NP e.g. at exchange:
      //NP      forces would be initialized with 0
      //NP      positions would be initialized from buffer data

      //NP restart always pulls from buffer
      if(op.operation() == OPERATION_RESTART)
          return true;

      //NP forward and reverse comm always pull from buffer
      //NP (thats why they are done)
      if(op.operation() == OPERATION_COMM_FORWARD ||
         op.operation() == OPERATION_COMM_REVERSE )
        return true;


      //NP exchange() and borders()
      if(op.operation() == OPERATION_COMM_BORDERS ||
         op.operation() == OPERATION_COMM_EXCHANGE )
      {
          //NP comm none and comm reverse dont pull from buffer
          if(communicationType_ == COMM_TYPE_NONE ||
             communicationType_ == COMM_TYPE_REVERSE)
             return false;

          //NP all others do
          return true;
      }

      // default
      return true;
  }

  /* ----------------------------------------------------------------------
   decide if unpack creates new element or overwrites existing data
  ------------------------------------------------------------------------- */

  inline bool ParticleDataContainerProperties::decideCreateNewElements(OperationProperties &op)
  {
      //NP have to decide at unpack if new elements are created or
      //NP existing ones are over-written

      //NP restart always creates new elements
      if(op.operation() == OPERATION_RESTART)
          return true;

      //NP exchange() and borders() always create new elements
      if(op.operation() == OPERATION_COMM_BORDERS ||
         op.operation() == OPERATION_COMM_EXCHANGE )
        return true;

      //NP forward and reverse comm never create new elements
      if(op.operation() == OPERATION_COMM_FORWARD ||
         op.operation() == OPERATION_COMM_REVERSE )
        return false;

      // default
      return false;
  }

  /* ----------------------------------------------------------------------
   fast test for reference frame
   note that rotation is only carried out for LEN_VEC==3
  ------------------------------------------------------------------------- */

    bool ParticleDataContainerProperties::isScaleInvariant()
    {
       return ( refFrame_ == REF_FRAME_INVARIANT ||
                refFrame_ == REF_FRAME_SCALE_TRANS_INVARIANT);
    }

    bool ParticleDataContainerProperties::isTranslationInvariant()
    {
        return ( refFrame_ == REF_FRAME_INVARIANT ||
                 refFrame_ == REF_FRAME_TRANS_ROT_INVARIANT ||
                 refFrame_ == REF_FRAME_SCALE_TRANS_INVARIANT ||
                 refFrame_ == REF_FRAME_TRANS_INVARIANT);
    }

    bool ParticleDataContainerProperties::isRotationInvariant()
    {
        return ( refFrame_ == REF_FRAME_INVARIANT ||
                 refFrame_ == REF_FRAME_TRANS_ROT_INVARIANT);
    }

  /* ----------------------------------------------------------------------
   ID operations
  ------------------------------------------------------------------------- */
  /*
  inline void ParticleDataContainerProperties::id(char *_id)
  {
      strcpy(_id,id_);
  }*/

  inline bool ParticleDataContainerProperties::matches_id(const char *_id)
  {
      if(strcmp(_id,id_) == 0) return true;
      return false;
  }

#endif
