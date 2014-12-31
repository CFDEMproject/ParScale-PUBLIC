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

#include "particle_data_container_properties.h"

#include "stdio.h"

using namespace PASCAL_NS;

    ParticleDataContainerProperties::ParticleDataContainerProperties()
    : id_(0),
      communicationType_(COMM_TYPE_MANUAL),
      refFrame_(REF_FRAME_UNDEFINED),
      restartType_(RESTART_TYPE_UNDEFINED),
      couplingType_(COUPLING_TYPE_UNDEFINED),
      do_read_(true),
      do_output_(false),
      needsBcast_(false),
      element_property_(false),
      scope_(0),
      scalePower_(-1),
      err_(false)
    {
    }

    ParticleDataContainerProperties::ParticleDataContainerProperties(const char *_id)
    : id_(0),
      communicationType_(COMM_TYPE_MANUAL),
      refFrame_(REF_FRAME_UNDEFINED),
      restartType_(RESTART_TYPE_UNDEFINED),
      couplingType_(COUPLING_TYPE_UNDEFINED),
      do_read_(true),
      do_output_(false),
      needsBcast_(false),
      element_property_(false),
      scope_(0),
      scalePower_(-1),
      err_(false)
    {
              if(id_) delete []id_;
              id_ = new char[strlen(_id)+1];
              strcpy(id_,_id);
    }


    ParticleDataContainerProperties::ParticleDataContainerProperties(const char *_id, const char* _comm,
               const char* _ref, const char *_restart,const char* _coupling,
               const char* _do_read, const char* _do_output,const char *_element_property,
               const char *_scope,int _scalePower
                                         )
    : id_(0),
      communicationType_(COMM_TYPE_MANUAL),
      refFrame_(REF_FRAME_UNDEFINED),
      restartType_(RESTART_TYPE_UNDEFINED),
      couplingType_(COUPLING_TYPE_UNDEFINED),
      do_read_(true),
      do_output_(false),
      element_property_(false),
      scope_(0),
      scalePower_(-1),
      err_(false)
    {
              if(id_) delete []id_;
              id_ = new char[strlen(_id)+1];
              strcpy(id_,_id);

              if      (strcmp(_comm,"comm_forward") == 0) communicationType_ = COMM_TYPE_FORWARD;
              else if (strcmp(_comm,"comm_forward_from_frame") == 0) communicationType_ = COMM_TYPE_FORWARD_FROM_FRAME;
              else if (strcmp(_comm,"comm_reverse") == 0) communicationType_ = COMM_TYPE_REVERSE;
              else if (strcmp(_comm,"comm_exchange_borders") == 0) communicationType_ = COMM_EXCHANGE_BORDERS;
              else if (strcmp(_comm,"comm_none") == 0) communicationType_ = COMM_TYPE_NONE;
              else if (strcmp(_comm,"comm_manual") == 0) communicationType_ = COMM_TYPE_MANUAL;
              else {printf("A"); communicationType_ = COMM_TYPE_UNDEFINED;}

              if      (strcmp(_ref,"frame_invariant") == 0) refFrame_ = REF_FRAME_INVARIANT;
              else if (strcmp(_ref,"frame_trans_rot_invariant") == 0) refFrame_ = REF_FRAME_TRANS_ROT_INVARIANT;
              else if (strcmp(_ref,"frame_scale_trans_invariant") == 0) refFrame_ = REF_FRAME_SCALE_TRANS_INVARIANT;
              else if (strcmp(_ref,"frame_trans_invariant") == 0) refFrame_ = REF_FRAME_TRANS_INVARIANT;
              else if (strcmp(_ref,"frame_general") == 0) refFrame_ = REF_FRAME_GENERAL;
              else {printf("B"); refFrame_ = REF_FRAME_UNDEFINED;}

              if      (strcmp(_restart,"restart_yes") == 0) restartType_ = RESTART_TYPE_YES;
              else if (strcmp(_restart,"restart_no") == 0) restartType_ = RESTART_TYPE_NO;
              else {printf("C"); restartType_ = RESTART_TYPE_UNDEFINED;}

              if      (strcmp(_coupling,"coupling_pull") == 0) couplingType_ = COUPLING_TYPE_PULL;
              else if (strcmp(_coupling,"coupling_push") == 0) couplingType_ = COUPLING_TYPE_PUSH;
              else if (strcmp(_coupling,"coupling_pull_push") == 0) couplingType_ = COUPLING_TYPE_PULL_PUSH;
              else if (strcmp(_coupling,"coupling_push_min_max") == 0) couplingType_ = COUPLING_TYPE_PUSH_MIN_MAX;
              else {printf("WARNING: you are using an undefined coupling type \n"); couplingType_ = COUPLING_TYPE_UNDEFINED;}

              if      (strcmp(_do_read,"read_yes") == 0) do_read_ = true;
              else if (strcmp(_do_read,"read_no") == 0) do_read_ = false;
              else err_ = true;

              if      (strcmp(_do_output,"output_yes") == 0) do_output_ = true;
              else if (strcmp(_do_output,"output_no") == 0) do_output_ = false;
              else err_ = true;

              if      (strcmp(_element_property,"element_property") == 0) element_property_ = true;
              else if (strcmp(_element_property,"global_property") == 0)  element_property_ = false;
              else err_ = true;

              if(_scope)
              {
                // filename might be name of the model
                int n = strlen(_scope) + 1;
                scope_ = new char[n];
                strcpy(scope_,_scope);
              }
              else
              {
                // default filename is name of the property
                int n = strlen(_id) + 1;
                scope_ = new char[n];
                strcpy(scope_,_id);
              }

              scalePower_ = _scalePower;
    }

    /* ----------------------------------------------------------------------
       destructor
    ------------------------------------------------------------------------- */

    ParticleDataContainerProperties::~ParticleDataContainerProperties()
    {
        delete [] scope_;
    }


    /* ----------------------------------------------------------------------
       set comm and reference properties
    ------------------------------------------------------------------------- */

    void ParticleDataContainerProperties::setProperties(ParticleDataContainerProperties &cp)
    {
        if(id_) delete []id_;
        id_ = new char[strlen(cp.id_)+1];
        strcpy(id_,cp.id_);

        if(scope_) delete []scope_;
        scope_ = new char[strlen(cp.scope_)+1];
        strcpy(scope_,cp.scope_);

        communicationType_  = cp.communicationType_;
        refFrame_           = cp.refFrame_;
        restartType_        = cp.restartType_;
        couplingType_       = cp.couplingType_;
        do_read_            = cp.do_read_;
        do_output_          = cp.do_output_;
        element_property_   = cp.element_property_;
        scalePower_         = cp.scalePower_;
    }

    bool ParticleDataContainerProperties::propertiesSetCorrectly()
    {
      if(refFrame_ == REF_FRAME_UNDEFINED ||
        communicationType_ == COMM_TYPE_UNDEFINED ||
        restartType_ == RESTART_TYPE_UNDEFINED ||
        couplingType_ == COUPLING_TYPE_UNDEFINED ||
        scalePower_ < 0)
            return false;

      return true;
    }
