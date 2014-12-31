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

#ifndef PSC_PARTICLE_DATA_CONTAINER_PROPERTIES_H
#define PSC_PARTICLE_DATA_CONTAINER_PROPERTIES_H

#include "string.h"
#include "stdio.h"

namespace PASCAL_NS
{

  /* ----------------------------------------------------------------------
   definition of reference frames, comm types etc
  ------------------------------------------------------------------------- */

  // buffer operation types (for push and pop)
  enum{ OPERATION_COMM_EXCHANGE,
        OPERATION_COMM_BCAST,
        OPERATION_COMM_BORDERS,
        OPERATION_COMM_FORWARD,
        OPERATION_COMM_REVERSE,
        OPERATION_RESTART,
        OPERATION_READ, // from JSON files
        OPERATION_OUTPUT,
        OPERATION_PULL,
        OPERATION_PUSH,
        OPERATION_UNDEFINED};

  // reference frame types
  // invariant: invariant to scaling, translation, rotation
  // trans invariant: invariant to translation, not invariant to scaling, rotation
  // trans+rot invariant: invariant to translation, rotation, not invariant to scaling
  // general: not invariant to scaling, translation, rotation

  enum{ REF_FRAME_UNDEFINED,
        REF_FRAME_INVARIANT,
        REF_FRAME_SCALE_TRANS_INVARIANT,
        REF_FRAME_TRANS_ROT_INVARIANT,
        REF_FRAME_TRANS_INVARIANT,
        REF_FRAME_GENERAL};

  // communication types

  enum{ // communication invoked manually
        COMM_TYPE_MANUAL,
        // only exchange and borders comm
        COMM_EXCHANGE_BORDERS,
        // forward comm every step
        COMM_TYPE_FORWARD,
        // forward comm based on reference frame setting
        // ie if mesh rotates, egdeVecs are communicated
        //NP does exchange, borders with buffer-initialized values
        COMM_TYPE_FORWARD_FROM_FRAME,
        // reverse comm every step
        //NP does exchange and borders with 0-initialized values
        COMM_TYPE_REVERSE,
        // no comm at all
        //NP does exchange and borders with 0-initialized values
        COMM_TYPE_NONE,
        // undefined state for error check
        COMM_TYPE_UNDEFINED};  // communication types

  // restart types

  enum{ RESTART_TYPE_UNDEFINED,
        RESTART_TYPE_YES,
        RESTART_TYPE_NO};

  // coupling types

  enum{ COUPLING_TYPE_UNDEFINED,
        COUPLING_TYPE_NONE,
        COUPLING_TYPE_PULL,
        COUPLING_TYPE_PUSH,
        COUPLING_TYPE_PUSH_MIN_MAX, //just pushes the first and last value in the array
        COUPLING_TYPE_PULL_PUSH};

  class OperationProperties
  {
      public:

        OperationProperties()
        : operation_(OPERATION_UNDEFINED),
          scale_(false),
          translate_(false),
          rotate_(false)
        {
        }
        OperationProperties(int _operation, bool _scale,bool _translate,bool _rotate)
        : operation_(_operation),
          scale_( _scale),
          translate_( _translate),
          rotate_( _rotate)
        {
        }

        inline void set_operation(int _operation)
        { operation_ = _operation; }

        inline int operation()
        { return operation_; }

        inline bool scale()
        { return scale_; }

        inline bool translate()
        { return translate_; }

        inline bool rotate()
        { return rotate_; }

      private:
        int operation_;
        bool scale_;
        bool translate_;
        bool rotate_;
  };


  class ParticleDataContainerProperties
  {
      public:

          ParticleDataContainerProperties();
          ParticleDataContainerProperties(const char *_id);
          ParticleDataContainerProperties(const char *_id, const char* _comm,
               const char* _ref, const char *_restart,const char* _coupling,
               const char* _do_read, const char* _do_output,const char *_element_property,
               const char *_scope,int _scalePower
                                         );
          ~ParticleDataContainerProperties();

          void setProperties(ParticleDataContainerProperties &cp);
          bool propertiesSetCorrectly();

          void needBCast()
          { needsBcast_ = true; }
          void needNoBCast()
          { needsBcast_ = false; }

          inline const char* id() const
          { return id_; }
          inline bool element_property() const
          { return element_property_; }
          inline const char* scope() const
          { return scope_; }
          inline bool matches_id(const char *_id);
          inline bool isScaleInvariant();
          inline bool isTranslationInvariant();
          inline bool isRotationInvariant();
          inline bool needsPush() const
          { 
            if ( 
                    COUPLING_TYPE_PUSH         == couplingType_ 
                 || COUPLING_TYPE_PULL_PUSH    == couplingType_ 
                 || COUPLING_TYPE_PUSH_MIN_MAX == couplingType_ 
               )
                return true;
            return false;
          }

          inline bool pushMax() const
          { 
            if (COUPLING_TYPE_PUSH_MIN_MAX == couplingType_ )
                return true;
            return false;
          }

          inline bool pushMin() const
          { 
            if (COUPLING_TYPE_PUSH_MIN_MAX == couplingType_ )
                return true;
            return false;
          }
          
          inline bool needsPull() const
          { 
            if (COUPLING_TYPE_PULL==couplingType_ || COUPLING_TYPE_PULL_PUSH==couplingType_)
                return true;
            return false;
          }


          inline bool needsRead() const
          { 
            return do_read_;
          }
    
          inline bool needsOutput() const
          { 
            return do_output_;
          }
          
          inline int scalePower()
          {  return scalePower_; }

          //NP decide on wheater at all an operation is performed here
          inline bool decidePackUnpackOperation(OperationProperties &op);

          //NP decide if operation performs data communication
          inline bool decideCommOperation(OperationProperties &op);

          //NP decide if unpack creates new element or overwrites existing data
          inline bool decideCreateNewElements(OperationProperties &op);

      private:

        char *id_;
        int  communicationType_;
        int  refFrame_;
        int  restartType_;
        int  couplingType_;
        bool do_read_;
        bool do_output_;
        bool needsBcast_;       // true if container needs bcast before using (eg because read from file)
        bool element_property_; // element property if true; global property if false
        char *scope_;           //JSON file where property is located
        int  scalePower_;
        bool err_;
  };

  // *************************************
  #include "particle_data_container_properties_I.h"
  // *************************************
}

#endif
