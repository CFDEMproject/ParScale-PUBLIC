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

#include "coupling.h"
#include "input.h"
#include "comm.h"
#include "output.h"
#include "stdlib.h"
#include "particle_data.h"
#include "particle_mesh.h"
#include "container_scalar.h"
#include "style_coupling_model.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   ModelBase Constructor & Destructor
------------------------------------------------------------------------- */

Coupling::Coupling(ParScale *ptr) :
    ParScaleBaseAccessible(ptr),
    couplingModel_map_(new std::map<std::string,CouplingModelCreator>()),
    couplingModels_(0)
{
    isInitialized_ = false;
    verbose_ = false;        
    // fill map with all models listed in style_model.h

    #define COUPLING_MODEL_CLASS
    #define CouplingModelStyle(key,Class) \
      (*couplingModel_map_)[#key] = &couplingModel_creator<Class>;
    #include "style_coupling_model.h"
    #undef CouplingModelStyle
    #undef COUPLING_MODEL_CLASS

    if(couplingModels_.size() > 0)
        error().throw_error_all(FLERR,"Currently supporting only one couplingmodel at a time\n");
}

Coupling::~Coupling()
{
    delete couplingModel_map_;
    couplingModels_.clear();
}

/* ----------------------------------------------------------------------
   one instance per model in style_model.h
------------------------------------------------------------------------- */
template <typename T>
CouplingModel *Coupling::couplingModel_creator(ParScale *ptr, char *name)
{
    return new T(ptr,name);
}

/* ----------------------------------------------------------------------
    Pull data from LIGGGHTS; gets called via ParScaleBaseInterfaceVector::init()
------------------------------------------------------------------------- */

void Coupling::read()
{
    for(unsigned int iModel=0; iModel < couplingModels_.size(); iModel++)
        couplingModels_[iModel].read();

    // pull all particle data
    // has to be done here so all particle data is available before
    // because has been bcasted / parallelized
    // data is either read or pulled
    particleData().pull();
}

/* ----------------------------------------------------------------------
    other stuff
------------------------------------------------------------------------- */

void Coupling::write()
{
}

void Coupling::bcast()
{
}

void Coupling::parallelize()
{
}

///////////////////////////////////////////////////////

void Coupling::init()
{

    for(int iModel=0; iModel < couplingModels_.size(); iModel++)
    {
        couplingModels_[iModel].init();
        if(couplingModels_[iModel].verbose())
            verbose_ = true;    

    }

    // error check
    external_code_in_control();

    isInitialized_ = true;
}

/* ----------------------------------------------------------------------
    re-arrange data,
    called from Driver::run() and ParScaleBaseInterfaceVector::init()
    called after Comm::pull()
    actual pull is done from particleData::pull()
------------------------------------------------------------------------- */
void Coupling::pull()
{
    // pull all particle data
    if(coupling().external_code_in_control())
    {
        // get number of local and global particles and re-alocate containers
        CustomValueTracker &cv = particleData().data();
        ContainerScalar<int> *id_cont = cv.getElementProperty<ContainerScalar<int> >("id");
        if(!id_cont)
            (ParScaleBase::error()).throw_error_all(FLERR,"Could not get element property 'id' needed by coupling.");

        int id_length_ = id_cont->size(); //the original length of the container. Important to detect change
        int id_lengthDelta_ = 0;


        couplingModel().pull_n_bodies(nlocal,nghost, nbody_all);

        double * tempIntraData_;
        if(verbose_)
        {
            //Allocate mem since needed for print out 
            tempIntraData_   = create<double>(tempIntraData_, particleMesh().nGridPoints());
            printf("\n[%d/%d]:***Coupling::pull()***\n", comm().me(),comm().nprocs());
            printf("[%d/%d]:current atom count: %d (local), %d (global).\n", 
                   comm().me(),comm().nprocs(),
                   cv.nbody(), cv.nbody_all());
            printf("[%d/%d]:atom count requested by coupling: %d (local), %d (ghost), %d (global).\n", 
                    comm().me(),comm().nprocs(),
                   nlocal, nghost, nbody_all);

            printf("[%d/%d]:received access to container with length %d and ids \n", 
                    comm().me(),comm().nprocs(),
                    id_cont->size()
                  );
            for(int iCont=0; iCont < (id_cont->size()); iCont++)
            {
                printf("[%d/%d]:id_cont[%d]: %d\n", 
                comm().me(),comm().nprocs(),
                iCont, id_cont->begin()[iCont]);
            }
        }

        if( cv.nbody_all() != nbody_all )
        {
            printf("[%d/%d]: cv.nbody_all(): %d,  nbody_all: %d \n", 
                    comm().me(),comm().nprocs(),
                    cv.nbody_all(),nbody_all
                  );
            (ParScaleBase::error()).throw_error_all(FLERR,"ParScale.nbody_all() != nbody_all from coupling! Check your coupling model, since the number of bodies in ParScale does not match that of the calling program.");
        }

        cv.grow_nbody(nlocal,nbody_all);
        id_lengthDelta_ = nlocal - id_length_; //get number of asynchronous atoms. positive if calling program has more
        id_length_      = id_cont->size();
        int *_id;
        _id = create<int>(_id,id_length_);
        memcpy(_id,id_cont->begin(),id_length_*sizeof(int));
        

        int ext_map_length     = 0;
        int *_ext_map = couplingModel().get_external_map(ext_map_length);

        if(_ext_map==NULL)
            (ParScaleBase::error()).throw_error_all(FLERR,"pointer to _ext_map not set!");

        //check current against external ids
        if(id_lengthDelta_>0)
        {
            int idsToBeDetected = id_lengthDelta_;
            int iExt=0;
            while( idsToBeDetected>0 && (iExt<ext_map_length) )
            {
                if( _ext_map[iExt] < nlocal ) //check if local ID is in range
                {
                    int extGlobalID = iExt + 1; //the global ID from the external map
                    bool foundGlobal=false;
                    for( int jLocal=0; jLocal<(nlocal-id_lengthDelta_); jLocal++ )//scan all old global ids for external global id
                    {
                        if( extGlobalID == _id[jLocal] )
                        {
                            foundGlobal = true;
                            continue;
                        }
                    }

                    if(!foundGlobal) //did not find global id, so need to add at the end
                    {
                        _id[nlocal-1+idsToBeDetected] = extGlobalID; //inser global Id at the end
                        idsToBeDetected--;
                    }
                    printf("[%d/%d]:couplingModel added extGlobalID %d into slot %d \n", 
                           comm().me(),comm().nprocs(),
                           extGlobalID, nlocal-1+idsToBeDetected
                           );
                }
                iExt++;
            }
        }
        

        if(verbose_)
        {
            printf("[%d/%d]:couplingModel attempts to pull %d (local), %d (ghost) and %d (global) bodies \n", 
                    comm().me(),comm().nprocs(),
                    nlocal, nghost, nbody_all
                  );
            printf("[%d/%d]:received external map of length %d with values: \n", 
                   comm().me(),comm().nprocs(),
                   ext_map_length
                  );
            for(int iGlobal=0; iGlobal < ext_map_length; iGlobal++)
            {
                printf("[%d/%d]:_ext_map[%d]: %d\n", 
                        comm().me(),comm().nprocs(),
                        iGlobal, _ext_map[iGlobal]
                      );
            }
            
            printf("[%d/%d]:created copy of id_container with length: %d\n", 
                    comm().me(),comm().nprocs(),
                    id_length_
                  );
        }

        if(verbose_)
        for(int idLocal=0; idLocal < id_length_; idLocal++)
        {
            particleData().setParticleIDPointer(0,idLocal);	
            particleData().returnIntraData(tempIntraData_);

            printf("[%d/%d]: BEFORESORTING idLocal: %d/(nlocal: %d, nghosts: %d), globalID: %d, intraData[0]: %g \n", 
                    comm().me(),comm().nprocs(),
                    idLocal, nlocal, nghost, _id[idLocal],tempIntraData_[0]
                  );
        }

        // use to sort data and afterwards destroy copy
        cv.sortPropsByExtMap(
                             _id, nlocal, 
                             id_length_,
                             _ext_map, ext_map_length,
                             verbose_, comm().me()
                            );

        //Check IDs. 
        //Must ensure that (i) global ids are unique for nlocal
        if(verbose_)
        for(int idLocal=0; idLocal < id_length_; idLocal++)
        {
            particleData().setParticleIDPointer(0,idLocal);	
            particleData().returnIntraData(tempIntraData_);

            printf("[%d/%d]: AFTERSORTING idLocal: %d/(nlocal: %d, nghosts: %d), globalID: %d, intraData[0]: %g  \n", 
                    comm().me(),comm().nprocs(),
                    idLocal, nlocal, nghost, _id[idLocal],tempIntraData_[0]
                  );

            delete [] tempIntraData_;
        }

        //Check external map
        for(int iLocal=0; iLocal < nlocal; iLocal++)
            if( _ext_map[_id[iLocal]-1]<0 )
                (ParScaleBase::error()).throw_error_all(FLERR,"List of globals ids invalid! _ext_map < 0");
            else if(  _ext_map[_id[iLocal]-1] >= nlocal)
                (ParScaleBase::error()).throw_error_all(FLERR,"List of globals ids invalid! _ext_map >= nlocal");
        
        //Finally pull the data        
        particleData().pull();

        destroy(_id);
        
    }

    if (strcmp(coupling().couplingModel().name(),"json") == 0)
    {
        particleData().pull();
    }

    couplingModels_[0].setPushPull();
}

/* ----------------------------------------------------------------------
    Push all particle data
    called from Driver::run()
------------------------------------------------------------------------- */

void Coupling::push()
{

    if(coupling().external_code_in_control())
    {
        if(verbose_)
            printf("\n***Coupling::push()***\n");
       
        //Since particles should be sorted, simply push all particle data
        particleData().push();
    }
    couplingModels_[0].setPushPull();
}

/* ----------------------------------------------------------------------
    Check if there is coupling to external code which is in control,
    i.e. the external code is the driver
------------------------------------------------------------------------- */

bool Coupling::external_code_in_control() const
{
    int ncontrol = 0;
    for(int iModel=0; iModel < couplingModels_.size(); iModel++)
    {
        if(couplingModels_[iModel].external_code_in_control())
            ncontrol++;
    }

    if(ncontrol > 1)
        (ParScaleBase::error()).throw_error_all(FLERR,"Can not use more than one coupling "
                                      "that controls ParScale");
    if(1 == ncontrol)
        return true;

    return false;
}

/* ----------------------------------------------------------------------
    Parse commands for Coupling class
------------------------------------------------------------------------- */

void Coupling::parse_command(int narg,char const* const* arg)
{
    //printf("narg %d\n",narg);
    int n = strlen(arg[0]) + 1;
    char *couplingModelName = new char[n];
    strcpy(couplingModelName,arg[0]);

    if (couplingModel_map_->find(arg[0]) != couplingModel_map_->end())
    {
        CouplingModelCreator model_creator = (*couplingModel_map_)[arg[0]];
        couplingModels_.push_back(model_creator(pascal_ptr(), couplingModelName));


        printf("...coupling model %s (ID: %lu) is registered with name %s\n\n",
               arg[0],
               couplingModels_.size()-1,
               couplingModels_[couplingModels_.size()-1].name());

    }
    else
        error().throw_error_all(FLERR,"Coupling PARSING: coupling model name not found: ",arg[0],"\n");

    delete [] couplingModelName;
}

/* ----------------------------------------------------------------------
    Functions that actually fill / dump the container data
    called for each container via ParticleData::pull() / push()

    currenlty only one coupling model active: EITHER JSON or LIGGGHTS
    TODO: extend in a way that the container knows where to fill itself from
------------------------------------------------------------------------- */

bool Coupling::fill_container_from_coupling(class ContainerBase &container) const
{
    // currently restriction to one coupling model
    return couplingModels_[0].fill_container_from_coupling(container);
}

/* ------------------------------------------------------------------------ */

bool Coupling::dump_container_to_coupling(class ContainerBase &container) const
{
 
   
    // currently restriction to one coupling model
    return couplingModels_[0].dump_container_to_coupling(container);
}
