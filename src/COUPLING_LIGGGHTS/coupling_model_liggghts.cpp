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

#include "coupling_model_liggghts.h"
#include "output.h"
#include "error.h"
#include "comm.h"
#include "container.h"
#include "coupling.h"
#include "particle_data.h"
#include "library.h" //include LIGGGHTS' C-style interface functions

#define BIGNEGNUMBER -1000000.

using namespace PASCAL_NS;
using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   CouplingModelLiggghts Constructor
------------------------------------------------------------------------- */

CouplingModelLiggghts::CouplingModelLiggghts(ParScale *ptr, const char *_name)
:
CouplingModel(ptr,_name)
{
    isInitialized_ = false;

    if(verbose_)
        output().write_screen_all("\n****** CouplingModelLiggghts: SETTING LIGGGHTS POINTER *********");
        
    lmp_ = static_cast<LAMMPS_NS::LAMMPS*>(callingProgram());
    if(lmp_==NULL)
        error().throw_error_all(FLERR,"CouplingModelLiggghts:: pointer to LIGGGHTS/LAMMPS object not set. Will abort.");
    else
      if(verbose_)
        output().write_screen_all("CouplingModelLiggghts:: pointer to LIGGGHTS/LAMMPS SET!");

    dumpName    = new char[30];
    type        = new char[30];
    pullName    = new char[90];
    pullName2   = new char[90];
}

/* ----------------------------------------------------------------------
   CouplingModelLiggghts Destructor
------------------------------------------------------------------------- */
CouplingModelLiggghts::~CouplingModelLiggghts()
{
    delete [] dumpName;
    delete [] type;
    delete [] pullName;
    delete [] pullName2;

}

/* ----------------------------------------------------------------------
   CouplingModelLiggghts Constructor
------------------------------------------------------------------------- */
void CouplingModelLiggghts::init()
{
    fix_coupling_ = static_cast<LAMMPS_NS::FixParScaleCouple*>(lmp_->modify->find_fix_style_strict("couple/pascal",0));
    if(fix_coupling_)
        output().write_screen_all("CouplingModelLiggghts has found fix of type 'couple/pascal'. Calling fix-functions now...");
    else
        error().throw_error_all(FLERR,"Pointer to LIGGGHTS fix 'couple/pascal' could not be set. Will abort.");

    exchangeEventsLocalId_          = fix_coupling_->exchangeEventsLocalId;
    exchangeEventsReceivingProcess_ = fix_coupling_->exchangeEventsReceivingProcess;

    isInitialized_ = true;
}

/* ----------------------------------------------------------------------
functions to get LIGGGHTS data
------------------------------------------------------------------------- */
void CouplingModelLiggghts::pull_n_bodies(int &_nlocal, int &_nghost, int &_nbody_all)
{
    //Pull-out values via LMP C-style library functions
    _nlocal     = *((int *) lammps_extract_global(lmp_,"nlocal"));
    _nbody_all = *((int *) lammps_extract_global(lmp_,"natoms"));
    _nghost    = *((int *) lammps_extract_global(lmp_,"nghost"));
    
    if(verbose_)
    {
        char msgstr[500];
        sprintf(msgstr, "CouplingModelLiggghts::pull_n_bodies, nlocal %d, nghost %d, natoms %d ",
                _nlocal, _nghost, _nbody_all);
        output().write_screen_all(msgstr);
    }
}

// - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - - 
void CouplingModelLiggghts::pull_box(double *_boxlo,double *_boxhi, double *_sublo,double *_subhi)
{
    //Copy pointers via LMP C-style library functions, just pull x since y,z are next in LIGGGHTS' mem
    memcpy(_boxlo,
           (double *)lammps_extract_global(lmp_,"boxxlo"), //pointer to x-value = 1st value in array,
           3*sizeof(double)
          );
    memcpy(_boxhi,
           (double *)lammps_extract_global(lmp_,"boxxhi"), //pointer to x-value = 1st value in array,
           3*sizeof(double)
          );

    memcpy(_sublo,
           (double *)lammps_extract_global(lmp_,"subxlo"), //pointer to x-value = 1st value in array,
           3*sizeof(double)
          );
    memcpy(_subhi,
           (double *)lammps_extract_global(lmp_,"subxhi"), //pointer to x-value = 1st value in array,
           3*sizeof(double)
          );
}

// - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - - 
void CouplingModelLiggghts::pull_proc_info(int *_procgrid,int *_myloc, int (&_procneigh)[3][2])
{

    //Copy pointers via LMP C-style library functions
    memcpy(_procgrid,
           (int *)lammps_extract_global(lmp_,"procx"),
           3*sizeof(int)
          );
    memcpy(_myloc,
           (int *)lammps_extract_global(lmp_,"mylocx"),
           3*sizeof(int)
          );

    //Get bare values
    _procneigh[0][0] = *((int *) lammps_extract_global(lmp_,"procneighxleft"));
    _procneigh[0][1] = *((int *) lammps_extract_global(lmp_,"procneighxright"));
    _procneigh[1][0] = *((int *) lammps_extract_global(lmp_,"procneighyleft"));
    _procneigh[1][1] = *((int *) lammps_extract_global(lmp_,"procneighyright"));
    _procneigh[2][0] = *((int *) lammps_extract_global(lmp_,"procneighzleft"));
    _procneigh[2][1] = *((int *) lammps_extract_global(lmp_,"procneighzright"));
}

// - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - - 
void CouplingModelLiggghts::pull_timeStepping_info(double &_deltaT, int &_neighAgo)
{
    //Pull-out values via LMP C-style library functions
    _deltaT     = *((double *) lammps_extract_global(lmp_,"dt"));
    memcpy(&_neighAgo,
           (int *)lammps_extract_global(lmp_,"ago"),
           sizeof(int)
          ); //use memcopy in order to have access to current value in caller program!
    
    if(verbose_)
    {
        char msgstr[500];
        sprintf(msgstr, "CouplingModelLiggghts::pull_timeStepping_info, deltaT %g, neighAgo %d.",
                _deltaT, _neighAgo);
        output().write_screen_all(msgstr);
    }

}
// - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - - 

int* CouplingModelLiggghts::get_external_map(int &length)
{
    if(fix_coupling_==NULL)
        error().throw_error_one(FLERR,"Pointer to LIGGGHTS fix 'couple/pascal' not set. Will abort.");

    return fix_coupling_->get_liggghts_map(length);
    
}

/* ----------------------------------------------------------------------
   pull/push for coupling model LIGGGHTS
------------------------------------------------------------------------- */
bool CouplingModelLiggghts::fill_container_from_coupling(class ContainerBase &container) const
{
 

    const ParticleDataContainerProperties& containerProps = container.prop();

    if(verbose_)
     printf("[%d/%d]:...fill container with id '%s', needsPull: %d, needsRead: %d double: %d, int: %d. \n",
               comm().me(), comm().nprocs(),
               containerProps.id(),
               containerProps.needsPull(),
               containerProps.needsRead(),
               container.isDoubleData(), container.isIntData()
            );
 
    if(!isInitialized_ || !containerProps.needsPull()) //need check here since read can trigger fill
        return false;

    if( particleData().nbody() > ( coupling().nlocal + coupling().nghost) )
        error().throw_error_one(FLERR,"Attempting to pull more bodies as available in the coupled code. This might cause seg faults.");

    if( containerProps.pullReset() && container.isIntData() )
    {
        printf( "***ERROR: vector property '%s' causes a problem*** \n", containerProps.id() );
        error().throw_error_one(FLERR,"CANNOT pullReset an integer data array. Not implemented/meaningful.");
    }

    if(verbose_)
    {
        output().write_screen_one("...container is now filled from LIGGGHTS");

        if(containerProps.pullReset())
            output().write_screen_one("...but only reset values will be pulled.");
    }


     //Check if pointer to Liggghts and the fix is present
    if(lmp_==NULL)
        error().throw_error_one(FLERR,"Pointer to LIGGGHTS/LAMMPS object not set. Will abort.");
    if(fix_coupling_==NULL)
        error().throw_error_one(FLERR,"Pointer to LIGGGHTS fix 'couple/pascal' not set. Will abort.");

    int len1,len2;
    if(1 == container.lenVec())
    {
        if( containerProps.pullReset() )
            error().throw_error_one(FLERR,"Cannot pullReset for a scalar container. Simply use pull.");
    
        if(verbose_)
        printf("[%d/%d]:...attempting to fill a scalar container with id '%s', double: %d, int: %d. \n",
               comm().me(), comm().nprocs(),
               containerProps.id(), container.isDoubleData(), container.isIntData());
              
        sprintf(type,     "%s","scalar-atom");


        void * ptr = fix_coupling_->find_push_property(containerProps.id(),type,len1,len2);
        if(ptr==NULL)
        {
             printf("***ERROR: scalar property '%s' causes a problem*** \n", containerProps.id());
             error().throw_error_one(FLERR,"could not find LIGGGHTS scalar push property.");
        }

        if(container.isDoubleData())
            memcpy((double*)container.begin_slow_dirty(),
                   (double*)ptr,
                   particleData().nbody()*sizeof(double)
                  );

        if(container.isIntData())
            memcpy((int*)container.begin_slow_dirty(),
                   (int*)ptr,
                   particleData().nbody()*sizeof(int)
                  );

        if(verbose_)
        {
            printf("...success! Received the following data:\n");
            for(int itemI=0;itemI<particleData().nbody(); itemI++)
            {

                if(container.isIntData())
                    printf("[%d/%d]:data[%d]: %d \n", 
                            comm().me(), comm().nprocs(),
                            itemI, 
                            static_cast<int*>(container.begin_slow_dirty())[itemI]
                          );
                if(container.isDoubleData())
                    printf("[%d/%d]:data[%d]: %g \n", 
                            comm().me(), comm().nprocs(),
                            itemI, 
                            static_cast<double*>(container.begin_slow_dirty())[itemI]
                          );

            }
            printf("\n");
        }                  
        return true;
    }
    else if( (3 == container.lenVec()) || (container.lenVec()>1 && containerProps.pullReset()) ) 
    {

        if(verbose_)
        printf("[%d/%d]:...attempting to fill a vector container with id '%s', double: %d, int: %d. \n",
               comm().me(), comm().nprocs(),
               containerProps.id(), container.isDoubleData(), container.isIntData());
              
        // create custom strings
        if( containerProps.pullReset() )
        {
            sprintf(type,     "%s","scalar-atom");
            sprintf(pullName, "%sReset", containerProps.id());
            sprintf(pullName2,"%sFluid", containerProps.id());
        }
        else
        {
            sprintf(type,     "%s","vector-atom");   
            sprintf(pullName, "%s", containerProps.id());
        }


        if(verbose_)
        printf("[%d/%d]:...attempting to dump a min value of a scalar container with id '%s', double: %d, int: %d. \n ...pullName(s) are %s/%s. \n",
               comm().me(), comm().nprocs(),
               containerProps.id(), container.isDoubleData(), container.isIntData(),  pullName,  pullName2);

        void * ptr;
        void * ptr2;
        ptr = fix_coupling_->find_pull_property(pullName,type,len1,len2);        

        if( containerProps.pullReset() )
        {
            ptr2 = fix_coupling_->find_pull_property(pullName2,type,len1,len2);        

            if( verbose_ && ptr==NULL )
               printf( "***WARNING: property '%s' cannot pull RESET scalar value  (name: %s) for reset \n", 
                       containerProps.id(), 
                       pullName
                     );
            if( verbose_ && ptr2==NULL )
               printf( "***WARNING: property '%s' cannot pull FLUID scalar value  (name: %s) for reset \n", 
                       containerProps.id(), 
                       pullName2
                     );

            if ( ptr==NULL && ptr2==NULL ) //if cannot find both, abort
                return true;
        }
        else
        {
            if(ptr==NULL)
            {
                printf( "***ERROR: vector property '%s' causes a problem*** \n", containerProps.id() );
                error().throw_error_one(FLERR,"could not find LIGGGHTS vector pull property.");
            }
        }

        if( containerProps.pullReset() )
        {

            //This routine ALWAYS assumes that the container isDoubleData   
            if( ptr!=NULL )  //pull and set RESET values,
            {
                double* values; //ptr to intra-particle data values for each particle
                int maxDataLocation = container.lenVecUsed()-2;
                for(int itemI=0;itemI<particleData().nbody(); itemI++)
                {
                    if( ((double*)ptr)[itemI] > BIGNEGNUMBER )
                    {
                       values = &(((double***)container.begin_slow_dirty())[itemI][0][0]); //ptr to first entry
                       //loop array and reset all 
                       for(int jInt=0; jInt<=maxDataLocation; jInt++)
                          values[jInt] = ((double*)ptr)[itemI];
                    }
                }

            }
            if( ptr2!=NULL ) //pull and set FLUID values
            {
                double* values; //ptr2 to FLUID data values for each particle
                int maxDataLocation = container.lenVecUsed()-1;
                for(int itemI=0;itemI<particleData().nbody(); itemI++)
                {
                    if(  ((double*)ptr2)[itemI] > BIGNEGNUMBER )
                    {
                       values    = &(((double***)container.begin_slow_dirty())[itemI][0][maxDataLocation]); //ptr to last entry = FLUID entry
                       values[0] = ((double*)ptr2)[itemI];
                    }
                }
            }
            return true;
        }
        else
        { //begin section standard vector data
          if(container.isDoubleData())
            for(int itemI=0;itemI<particleData().nbody(); itemI++)
              memcpy(&(((double**)container.begin_slow_dirty())[itemI]),
                     &(((double**)ptr)[itemI]),
                     3*sizeof(double)
                    );

          if(container.isIntData())
            for(int itemI=0;itemI<particleData().nbody(); itemI++)
              memcpy(&(((int**)container.begin_slow_dirty())[itemI]),
                     &(((int**)ptr)[itemI]),
                     3*sizeof(int)
                    );

       
          if(verbose_)
          {
            printf("...success! Received the following data:\n");
            for(int itemI=0;itemI<particleData().nbody(); itemI++)
            {

                if(container.isIntData())
                    printf("[%d/%d]:int data[%d]: %d %d %d \n", 
                            comm().me(), comm().nprocs(),
                            itemI, 
                            static_cast<int**>(container.begin_slow_dirty())[itemI][0],
                            static_cast<int**>(container.begin_slow_dirty())[itemI][1],
                            static_cast<int**>(container.begin_slow_dirty())[itemI][2]
                          );
                          
                if(container.isDoubleData())
                    printf("[%d/%d]:double data[%d]: %g %g %g \n", 
                            comm().me(), comm().nprocs(),
                            itemI, 
                            static_cast<double**>(container.begin_slow_dirty())[itemI][0],
                            static_cast<double**>(container.begin_slow_dirty())[itemI][1],
                            static_cast<double**>(container.begin_slow_dirty())[itemI][2]
                          );

            }
            printf("\n");
          }    
                          
          return true;
        } //end section standard vector data
    }
    else
        error().throw_error_one(FLERR,"Container.lenVec not in range, or memcpy not implemented. Aborting...");
    if (1 < container.nVec())
        error().throw_error_one(FLERR,"1 < container.nVec(). Aborting...");
        

    // should be if(isFilled_ || !this->containerProperties_.decidePackUnpackOperation(op)) there
    // check if it actually works for (a) read from file, (b) coupling
    error().throw_error_one(FLERR,"fill_container_from_coupling finished ILLEGALY. Abort.");

    return false;
}

///////////////////////////////////////////////////////////////////////////////////////////7
bool CouplingModelLiggghts::dump_container_to_coupling(class ContainerBase &container) const
{

    const ParticleDataContainerProperties& containerProps = container.prop();

    if(verbose_)
     printf("[%d/%d]:...dump container with id '%s', needsPush: %d, double: %d, int: %d. \n",
               comm().me(), comm().nprocs(),
               containerProps.id(),
               containerProps.needsPush(),
               container.isDoubleData(), container.isIntData()
            );
 
    if(!isInitialized_ || !containerProps.needsPush())
        return false;

    if(verbose_)
    {
        output().write_screen_one("...container is now dumped to LIGGGHTS");
        if(containerProps.pushMax())
            output().write_screen_one("...but max value will be dumped.");
        if(containerProps.pushMin())
            output().write_screen_one("...but min value will be dumped.");
    }   

     //Check if pointer to Liggghts and the fix is present
    if(lmp_==NULL)
        error().throw_error_one(FLERR,"Pointer to LIGGGHTS/LAMMPS object not set. Will abort.");
    if(fix_coupling_==NULL)
        error().throw_error_one(FLERR,"Pointer to LIGGGHTS fix 'couple/pascal' not set. Will abort.");

    int len1,len2;
    if(containerProps.pushMin())
    {
        sprintf(type,     "%s","scalar-atom");

        void * ptr;
        int dataLocation = 0;
        if( strcmp(containerProps.id(),"heat") == 0)
              sprintf(dumpName,"%sMin", "Temp");
        else
              sprintf(dumpName,"%sMin", containerProps.id());

        if(verbose_)
        printf("[%d/%d]:...attempting to dump a min value of a scalar container with id '%s', double: %d, int: %d. \n ...dumpName is %s. \n",
               comm().me(), comm().nprocs(),
               containerProps.id(), container.isDoubleData(), container.isIntData(),  dumpName);
        ptr = fix_coupling_->find_pull_property(dumpName,type,len1,len2);     

        if(ptr==NULL)
            error().throw_error_one(FLERR,"could not find LIGGGHTS scalar pull property.");

        if(container.isDoubleData())
          for(int itemI=0;itemI<particleData().nbody(); itemI++)
            memcpy(&(((double*)ptr)[itemI]),
                   &(((double***)container.begin_slow_dirty())[itemI][0][dataLocation]),
                    sizeof(double)
                  );
        if(container.isIntData())
          for(int itemI=0;itemI<particleData().nbody(); itemI++)
            memcpy(&(((int*)ptr)[itemI]),
                   &(((int***)container.begin_slow_dirty())[itemI][0][dataLocation]),
                    sizeof(int)
                  );

        if(verbose_)
        {
            printf("...success! Pushed the following data:\n");
            for(int itemI=0;itemI<particleData().nbody(); itemI++)
            {

                if(container.isDoubleData())
                    printf("[%d/%d]:double data[%d]: %g \n", 
                            comm().me(), comm().nprocs(),
                            itemI, 
                            static_cast<double*>(ptr)[itemI]
                          );
                if(container.isIntData())
                    printf("[%d/%d]:int data[%d]: %d \n", 
                            comm().me(), comm().nprocs(),
                            itemI, 
                            static_cast<int*>(ptr)[itemI]
                          );

            }
            printf("\n");
        }
    }
    if(containerProps.pushMax())
    {

        sprintf(type,     "%s","scalar-atom");
        void * ptr;
        int dataLocation = container.lenVecUsed()-2;

        if( strcmp(containerProps.id(),"heat") == 0)
           sprintf(dumpName,"%s", "Temp");
        else
           sprintf(dumpName,"%sMax", containerProps.id());

        if(verbose_)
        printf("[%d/%d]:...attempting to dump a max value of a scalar container with id '%s', double: %d, int: %d. \n...dumpName is %s. \n",
               comm().me(), comm().nprocs(),
             containerProps.id(), container.isDoubleData(), container.isIntData(),  dumpName);
        ptr = fix_coupling_->find_pull_property(dumpName,type,len1,len2);   

        if(ptr==NULL)
            error().throw_error_one(FLERR,"could not find LIGGGHTS scalar pull property.");

        if(container.isDoubleData())
          for(int itemI=0;itemI<particleData().nbody(); itemI++)
            memcpy(&(((double*)ptr)[itemI]),
                   &(((double***)container.begin_slow_dirty())[itemI][0][dataLocation]),
                    sizeof(double)
                  );

        if(container.isIntData())
          for(int itemI=0;itemI<particleData().nbody(); itemI++)
            memcpy(&(((int*)ptr)[itemI]),
                   &(((int***)container.begin_slow_dirty())[itemI][0][dataLocation]),
                    sizeof(int)
                  );

        if(verbose_)
        {
            printf("...success! Pushed the following data:\n");
            for(int itemI=0;itemI<particleData().nbody(); itemI++)
            {

                if(container.isDoubleData())
                    printf("[%d/%d]:double data[%d]: %g \n", 
                            comm().me(), comm().nprocs(),
                            itemI, 
                            static_cast<double*>(ptr)[itemI]
                          );
                if(container.isIntData())
                    printf("[%d/%d]:int data[%d]: %d \n", 
                            comm().me(), comm().nprocs(),
                            itemI, 
                            static_cast<int*>(ptr)[itemI]
                          );

            }
            printf("\n");
        }
        return true;
    }
    else if(1 == container.lenVec())
    {
    
        if(verbose_)
        printf("[%d/%d]:...attempting to dump a scalar container with id '%s', double: %d, int: %d. \n",
               comm().me(), comm().nprocs(),
               containerProps.id(), container.isDoubleData(), container.isIntData());
              
        sprintf(type,     "%s","scalar-atom");
        void * ptr;

        ptr = fix_coupling_->find_pull_property(containerProps.id(),type,len1,len2);
        if(ptr==NULL)
        {   
            printf("[%d/%d]: problem with containerProps.id(): %s. \n",
               comm().me(), comm().nprocs(),
               containerProps.id());
            error().throw_error_one(FLERR,"could not find LIGGGHTS scalar pull property.");
        }

            
        if(container.isDoubleData())
            memcpy((double*)ptr,
                   (double*)container.begin_slow_dirty(),
                   particleData().nbody()*sizeof(double)
                  );

        if(container.isIntData())
            memcpy((int*)ptr,
                   (int*)container.begin_slow_dirty(),
                   particleData().nbody()*sizeof(int)
                  );

        if(verbose_)
        {
            printf("...success! Pushed the following data:\n");
            for(int itemI=0;itemI<particleData().nbody(); itemI++)
            {

                if(container.isDoubleData())
                    printf("[%d/%d]:double data[%d]: %g \n", 
                            comm().me(), comm().nprocs(),
                            itemI, 
                            static_cast<double*>(ptr)[itemI]
                          );
                if(container.isIntData())
                    printf("[%d/%d]:int data[%d]: %d \n", 
                            comm().me(), comm().nprocs(),
                            itemI, 
                            static_cast<int*>(ptr)[itemI]
                          );

            }
            printf("\n");
        }
        return true;
    }
    else if(3 == container.lenVec()) //is vector!
    {
    
        if(verbose_)
        printf("[%d/%d]:...attempting to dump a vector container with id '%s', double: %d, int: %d. \n",
              comm().me(), comm().nprocs(),
              containerProps.id(), container.isDoubleData(), container.isIntData());
              
        sprintf(type,     "%s","vector-atom");
        void * ptr = fix_coupling_->find_pull_property(containerProps.id(),type,len1,len2);
        if(ptr==NULL)
        {
            printf("***ERROR: scalar property '%s' causes a problem*** \n", containerProps.id());
            error().throw_error_one(FLERR,"could not find LIGGGHTS vector pull property.");
        }
        
        if(container.isDoubleData())
          for(int itemI=0;itemI<particleData().nbody(); itemI++)
            memcpy(&(((double**)ptr)[itemI]),
                   &(((double**)container.begin_slow_dirty())[itemI]),
                   3*sizeof(double)
                  );

        if(container.isIntData())
          for(int itemI=0;itemI<particleData().nbody(); itemI++)
            memcpy(&(((int**)ptr)[itemI]),
                   &(((int**)container.begin_slow_dirty())[itemI]),
                   3*sizeof(int)
                  );

        if(verbose_)
        {
            printf("...success! Pushed the following data:\n");
            for(int itemI=0;itemI<particleData().nbody(); itemI++)
            {

                if(container.isDoubleData())
                    printf("[%d/%d]:double data[%d]: %g %g %g \n", 
                            comm().me(), comm().nprocs(),
                            itemI, 
                            static_cast<double**>(ptr)[itemI][0],
                            static_cast<double**>(ptr)[itemI][1],
                            static_cast<double**>(ptr)[itemI][2]
                          );

                if(container.isIntData())
                    printf("[%d/%d]:int data[%d]: %d %d %d \n", 
                            comm().me(), comm().nprocs(),
                            itemI, 
                            static_cast<int**>(ptr)[itemI][0],
                            static_cast<int**>(ptr)[itemI][1],
                            static_cast<int**>(ptr)[itemI][2]
                          );
                          
            }
            printf("\n");
        }    
                          
        return true;
    }
    else
        error().throw_error_one(FLERR,"Container.lenVec not in range, or memcpy not implemented. Aborting...");
    if (1 < container.nVec())
        error().throw_error_one(FLERR,"1 < container.nVec(). Aborting...");
        

    error().throw_error_one(FLERR,"dump_container_to_coupling finished ILLEGALY. Abort.");
    
    return false;
}
