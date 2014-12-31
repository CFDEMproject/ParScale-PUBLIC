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
#include "container.h"
#include "particle_data.h"
#include "library.h" //include LIGGGHTS' C-style interface functions

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
        output().write_screen_one("\n****** CouplingModelLiggghts: SETTING LIGGGHTS POINTER *********");
        
    lmp_ = static_cast<LAMMPS_NS::LAMMPS*>(callingProgram());
    if(lmp_==NULL)
        error().throw_error_one(FLERR,"CouplingModelLiggghts:: pointer to LIGGGHTS/LAMMPS object not set. Will abort.");
    else
      if(verbose_)
        output().write_screen_one("CouplingModelLiggghts:: pointer to LIGGGHTS/LAMMPS SET!");

}

/* ----------------------------------------------------------------------
   CouplingModelLiggghts Constructor
------------------------------------------------------------------------- */
void CouplingModelLiggghts::init()
{
    fix_coupling_ = static_cast<LAMMPS_NS::FixParScaleCouple*>(lmp_->modify->find_fix_style_strict("couple/pascal",0));
    if(fix_coupling_)
        printf("CouplingModelLiggghts has found fix of type 'couple/pascal'. Calling fix-functions now... \n");
    else
        error().throw_error_one(FLERR,"Pointer to LIGGGHTS fix 'couple/pascal' could not be set. Will abort.");

    isInitialized_ = true;
}

/* ----------------------------------------------------------------------
functions to get LIGGGHTS data
------------------------------------------------------------------------- */
void CouplingModelLiggghts::pull_n_bodies(int &_nbody, int &_nbody_all)
{
    //Pull-out values via LMP C-style library functions
    _nbody     = *((int *) lammps_extract_global(lmp_,"nlocal"));
    _nbody_all = *((int *) lammps_extract_global(lmp_,"natoms"));
}

// - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - - 
void CouplingModelLiggghts::pull_box(double *_boxlo,double *_boxhi, double *_sublo,double *_subhi)
{
    //Get pointers via LMP C-style library functions
    _boxlo     = (double *) lammps_extract_global(lmp_,"boxxlo"); //pointer to x-value = 1st value in array
    _boxhi     = (double *) lammps_extract_global(lmp_,"boxxhi"); //pointer to x-value = 1st value in array
    _sublo     = (double *) lammps_extract_global(lmp_,"subxlo"); //pointer to x-value = 1st value in array
    _subhi     = (double *) lammps_extract_global(lmp_,"subxhi"); //pointer to x-value = 1st value in array
}

// - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - -  - - - - - - 
void CouplingModelLiggghts::pull_proc_info(int *_procgrid,int *_myloc, int (&_procneigh)[3][2])
{

    //Get pointers via LMP C-style library functions
    _procgrid  = (int *) lammps_extract_global(lmp_,"procx"); //pointer to x-value = 1st value in array
    _myloc     = (int *) lammps_extract_global(lmp_,"mylocx"); //pointer to x-value = 1st value in array

    _procneigh[0][0] = *((int *) lammps_extract_global(lmp_,"procneighxleft"));
    _procneigh[0][1] = *((int *) lammps_extract_global(lmp_,"procneighxright"));
    _procneigh[1][0] = *((int *) lammps_extract_global(lmp_,"procneighyleft"));
    _procneigh[1][1] = *((int *) lammps_extract_global(lmp_,"procneighyright"));
    _procneigh[2][0] = *((int *) lammps_extract_global(lmp_,"procneighzleft"));
    _procneigh[2][1] = *((int *) lammps_extract_global(lmp_,"procneighzright"));
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
     printf("...processing container with id '%s', needsPull: %d, needsRead: %d double: %d, int: %d. \n",
               containerProps.id(),
               containerProps.needsPull(),
               containerProps.needsRead(),
               container.isDoubleData(), container.isIntData()
            );
 
    if(!isInitialized_ || !containerProps.needsPull()) //need check here since read can trigger fill
        return false;

    if(verbose_)
        output().write_screen_one("...container is now filled from LIGGGHTS");

     //Check if pointer to Liggghts and the fix is present
    if(lmp_==NULL)
        error().throw_error_one(FLERR,"Pointer to LIGGGHTS/LAMMPS object not set. Will abort.");
    if(fix_coupling_==NULL)
        error().throw_error_one(FLERR,"Pointer to LIGGGHTS fix 'couple/pascal' not set. Will abort.");

    int len1,len2;
    if(1 == container.lenVec())
    {
    
        if(verbose_)
        printf("...attempting to fill a scalar container with id '%s', double: %d, int: %d. \n",
               containerProps.id(), container.isDoubleData(), container.isIntData());
              
        const char*type = "scalar-atom";
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
                    printf("data[%d]: %d \n", 
                            itemI, 
                            static_cast<int*>(container.begin_slow_dirty())[itemI]
                          );
                if(container.isDoubleData())
                    printf("data[%d]: %g \n", 
                            itemI, 
                            static_cast<double*>(container.begin_slow_dirty())[itemI]
                          );

            }
            printf("\n");
        }                  
        return true;
    }
    else if(3 == container.lenVec()) 
    {
        error().throw_error_one(FLERR,"Pulling of vector quantities not functional!");

        if(verbose_)
        printf("...attempting to fill a vector container with id '%s', double: %d, int: %d. \n",
              containerProps.id(), container.isDoubleData(), container.isIntData());
              
        const char*type = "vector-atom";
        void * ptr = fix_coupling_->find_push_property(containerProps.id(),type,len1,len2);
        if(ptr==NULL)
        {
            printf("***ERROR: scalar property '%s' causes a problem*** \n", containerProps.id());
            error().throw_error_one(FLERR,"could not find LIGGGHTS vector push property.");
        }

        if(container.isDoubleData())
            memcpy((double**)container.begin_slow_dirty(),
                   (double**)ptr,
                   3*particleData().nbody()*sizeof(double)
                  );

        if(container.isIntData())
            memcpy((int**)container.begin_slow_dirty(),
                   (int**)ptr,
                   3*particleData().nbody()*sizeof(int)
                  );
       
        if(verbose_)
        {
            printf("...success! Received the following data:\n");
            for(int itemI=0;itemI<particleData().nbody(); itemI++)
            {

                if(container.isIntData())
                    printf("int data[%d]: %d %d %d \n", 
                            itemI, 
                            static_cast<int**>(container.begin_slow_dirty())[itemI][0],
                            static_cast<int**>(container.begin_slow_dirty())[itemI][1],
                            static_cast<int**>(container.begin_slow_dirty())[itemI][2]
                          );
                          
                if(container.isDoubleData())
                    printf("double data[%d]: %g %g %g \n", 
                            itemI, 
                            static_cast<double**>(container.begin_slow_dirty())[itemI][0],
                            static_cast<double**>(container.begin_slow_dirty())[itemI][1],
                            static_cast<double**>(container.begin_slow_dirty())[itemI][2]
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
     printf("...processing container with id '%s', needsPush: %d, double: %d, int: %d. \n",
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
        const char*type = "scalar-atom";
        char * dumpName = new char[30];
        void * ptr;
        int dataLocation = 0;
        if( strcmp(containerProps.id(),"heat") == 0)
              dumpName = "TempMin";
        else
              sprintf(dumpName,"%sMin", containerProps.id());

        if(verbose_)
        printf("...attempting to dump a min value of a scalar container with id '%s', double: %d, int: %d. \n ...dumpName is %s. \n",
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
                    printf("double data[%d]: %g \n", 
                            itemI, 
                            static_cast<double*>(ptr)[itemI]
                          );
                if(container.isIntData())
                    printf("int data[%d]: %d \n", 
                            itemI, 
                            static_cast<int*>(ptr)[itemI]
                          );

            }
            printf("\n");
        }
    }
    if(containerProps.pushMax())
    {

        const char*type = "scalar-atom";
        char * dumpName = new char[30];
        void * ptr;
        int dataLocation = container.lenVecUsed()-2;

        if( strcmp(containerProps.id(),"heat") == 0)
           dumpName = "Temp";
        else
           sprintf(dumpName,"%sMax", containerProps.id());

        if(verbose_)
        printf("...attempting to dump a max value of a scalar container with id '%s', double: %d, int: %d. \n...dumpName is %s. \n",
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
                    printf("double data[%d]: %g \n", 
                            itemI, 
                            static_cast<double*>(ptr)[itemI]
                          );
                if(container.isIntData())
                    printf("int data[%d]: %d \n", 
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
        printf("...attempting to dump a scalar container with id '%s', double: %d, int: %d. \n",
               containerProps.id(), container.isDoubleData(), container.isIntData());
              
        const char*type = "scalar-atom";
        void * ptr;

        ptr = fix_coupling_->find_pull_property(containerProps.id(),type,len1,len2);
        if(ptr==NULL)
            error().throw_error_one(FLERR,"could not find LIGGGHTS scalar pull property.");

            
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
                    printf("double data[%d]: %g \n", 
                            itemI, 
                            static_cast<double*>(ptr)[itemI]
                          );
                if(container.isIntData())
                    printf("int data[%d]: %d \n", 
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
        printf("...attempting to dump a vector container with id '%s', double: %d, int: %d. \n",
              containerProps.id(), container.isDoubleData(), container.isIntData());
              
        const char*type = "vector-atom";
        void * ptr = fix_coupling_->find_pull_property(containerProps.id(),type,len1,len2);
        if(ptr==NULL)
        {
            printf("***ERROR: scalar property '%s' causes a problem*** \n", containerProps.id());
            error().throw_error_one(FLERR,"could not find LIGGGHTS vector pull property.");
        }
        
        if(container.isDoubleData())
            memcpy((double**)ptr,
                   (double**)container.begin_slow_dirty(),
                   3*particleData().nbody()*sizeof(double)
                  );
                  
        if(container.isIntData())
            memcpy((int**)ptr,
                   (int**)container.begin_slow_dirty(),
                   3*particleData().nbody()*sizeof(int)
                  );

        if(verbose_)
        {
            printf("...success! Pushed the following data:\n");
            for(int itemI=0;itemI<particleData().nbody(); itemI++)
            {

                if(container.isDoubleData())
                    printf("double data[%d]: %g %g %g \n", 
                            itemI, 
                            static_cast<double**>(ptr)[itemI][0],
                            static_cast<double**>(ptr)[itemI][1],
                            static_cast<double**>(ptr)[itemI][2]
                          );

                if(container.isIntData())
                    printf("int data[%d]: %d %d %d \n", 
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
