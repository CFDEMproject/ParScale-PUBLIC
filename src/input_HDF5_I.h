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

int Input::write_containersHDF5(ContainerBase &container) const
{

    if(!writeHDF5_)
        return 0;

    if(!container.prop().needsOutput())
        return 0;

//    output().write_screen_one("Input::write_containersHDF5 is writing...");

    const char *id = container.prop().id();
    const char *scope = container.prop().scope();


    // if scope equals ID, try to open JSON file here directly
    if(strcmp(id,scope) == 0)
    {

        if(container.nVec()>1)
            error().throw_error_one(FLERR,"Cannot write if container.nVec()>1.");

        //generate filenames
        char filepath[100],hd5file[200];
        sprintf(filepath,"./%s/%.6f",
                runDirectory_,
                control().simulationState().time()
               );

        sprintf(filepath,"./%s/%.6f/HDF5",
                runDirectory_,
                control().simulationState().time()
               );

        if( comm().is_proc_0() )
        {
            mkdir(filepath, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IWOTH);
            mkdir(filepath, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IWOTH);
            sprintf(hd5file,"./%s/%.6f/HDF5/%s.h5",
                    runDirectory_,
                    control().simulationState().time(),
                    scope
                   );
        }

        comm().wait();

        //Skip this process in case it is empty
        if(container.size()==0)
            return 0;

        if( !comm().is_proc_0() )
            sprintf(hd5file,"./%s/%.6f/HDF5/%s.%d.h5",
                    runDirectory_,
                    control().simulationState().time(),
                    scope, comm().me()
                   );

        H5std_string    FILE_NAME( hd5file );
        H5std_string    DATASET_NAME( scope );

        const int     NX = container.size();
        const int     NY = container.lenVecUsed();                    // dataset dimensions
        const int     RANK = 2;

        //printf("Container %s should contain %i particles and have the length of %i \n", scope, NX, NY);

        double H5data[NX][NY];
        int    H5dataInt[NX][NY];
        //Get access to container content
        void *ptr = container.begin_slow_dirty();   //ptr to scalar per-particle data

        for(int iPart=0;iPart < container.size();iPart++) //Must only loop over used IDs!
        {
            for(int jVec=0;jVec < container.lenVecUsed();jVec++)
            {
                if(container.lenVec()>1)
                {
                   if(container.isDoubleData())
                       H5data[iPart][jVec]= static_cast<double***>(ptr)[iPart][0][jVec];
                   else if(container.isIntData())
                        H5dataInt[iPart][jVec] = static_cast<int***>(ptr)[iPart][0][jVec];
                   else
                        error().throw_error_one(FLERR,"Data format of container unknown.\n");
                }
                else
                {
                   if(container.isDoubleData())
                       H5data[iPart][jVec]= static_cast<double*>(ptr)[iPart];
                   else if(container.isIntData())
                        H5dataInt[iPart][jVec] = static_cast<int*>(ptr)[iPart];
                   else
                        error().throw_error_one(FLERR,"Data format of container unknown.\n");
                }
            }
        }

       try
       {

          Exception::dontPrint(); //print errors afterwards

          H5File* file = new H5File( FILE_NAME, H5F_ACC_TRUNC );    //Create a new file using H5F_ACC_TRUNC

          hsize_t  dimsf[] = {static_cast<hsize_t>(NX), 
                              static_cast<hsize_t>(NY)};                 // dataset dimensions
          //printf("dimsf for scope %s = [%i][%i] \n",scope, dimsf[0],dimsf[1]);
          //dimsf[0] = NX;                        // number of particles
          //dimsf[1] = NY;                        // length of container e.g. number of grid points
          DataSpace dataspace( RANK, dimsf);   // make space


          if(container.isDoubleData())
          {
              IntType datatype( PredType::NATIVE_DOUBLE );         //Define datatype for the data
              datatype.setOrder( H5T_ORDER_LE );
              DataSet* dataset = new DataSet(file->createDataSet( DATASET_NAME, datatype, dataspace )); //Create a new dataset within the file using
              dataset->write( H5data, PredType::NATIVE_DOUBLE ); //Write the data to the dataset using default memory space, file
              delete dataset;
          }
          else if(container.isIntData())
          {
              IntType datatype( PredType::NATIVE_INT );         //Define datatype for the data
              datatype.setOrder( H5T_ORDER_LE );
              DataSet* dataset = new DataSet(file->createDataSet( DATASET_NAME, datatype, dataspace )); //Create a new dataset within the file using
              dataset->write( H5dataInt, PredType::NATIVE_INT ); //Write the data to the dataset using
              delete dataset;
          }
          else
              error().throw_error_one(FLERR,"Data format of container unknown.\n");

          //printf("Wrote HDF5 Data for scope %s! \n",scope);

          delete file;
       }  // end of try block

       // catch failure caused by the H5File operations
       catch( FileIException error )
       {
          error.printError();
          return -1;
       }

       // catch failure caused by the DataSet operations
       catch( DataSetIException error )
       {
          error.printError();
          return -1;
       }

       // catch failure caused by the DataSpace operations
       catch( DataSpaceIException error )
       {
          error.printError();
          return -1;
       }

       // catch failure caused by the DataSpace operations
       catch( DataTypeIException error )
       {
          error.printError();
          return -1;
       }
      // printf("\n");
       return 0;  // successfully terminated
    }
    return 0;
}
