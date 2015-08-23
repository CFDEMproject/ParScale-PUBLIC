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


#include "particle_mesh.h"
#include "particle_data.h"
#include "input.h"
#include "output.h"
#include "control.h"
#include "comm.h"
#include "error.h"
#include "stdio.h"
#include "input_properties.h"
#include <iostream>
#include <iomanip>
#include <sys/stat.h>
using std::cout;
using std::endl;

#include <cstring>

const int MAX_CHARS_PER_LINE = 1024;
const int MAX_TOKENS_PER_LINE = 50;


using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   Constructor / Destructor
------------------------------------------------------------------------- */

Input::Input(ParScale *ptr) : ParScaleBase(ptr),
    infile_(0),
    myName("test"),
    writePrecision_(5)
{

    runDirectory_ = new char[512];
    sprintf(runDirectory_, "%s", "./");

    printf("default run directory: %s \n", runDirectory_);

    infileSet_ = false;

    process_JSONInput();

}

// * * * * * * * * * * * * * * * * * * * * *
Input::~Input()
{
    //if (infile_) infile_->close();
    infile_ = 0;

    delete [] runDirectory_;
}

/* ----------------------------------------------------------------------
   set the input script & run directory name
------------------------------------------------------------------------- */
void Input::set_input_script(const char *filename)
{
    myFile.open(filename, std::ifstream::in);
    infile_ = &myFile; //set pointer to infile

    infileSet_ = true;
   
}

// * * * * * * * * * * * * * * * * * * * * *
void Input::set_runDirectory(const char *dir)
{
    sprintf(runDirectory_, "%s", dir);
    printf("\nParScale will run in the '%s' directory ...\n", runDirectory_);
}

/* ----------------------------------------------------------------------
   process an input command
------------------------------------------------------------------------- */
void Input::process_input(const char * command)
{
    {
//        output().write_screen_all(" ParScale is processing a single input line ...");

        // process the tokens
        char* commandChar = strdup(command);

        const char* token[MAX_TOKENS_PER_LINE] = {}; // initialize to 0
        token[0] = strtok(commandChar, " \t"); // first token

        int n = 0; // a for-loop index
        for (n = 1; n < MAX_TOKENS_PER_LINE; n++)
        {
            token[n] = strtok(0, " \t"); // subsequent tokens
            if (!token[n]) break; // no more tokens
        }
        std::string str(token[0], strlen(token[0]));
        pascalInterface().parse_command(str,n-1,&(token[1]));

        free(commandChar);

//        output().write_screen_all("**Processed single input line.");
    }
    
}

/* ----------------------------------------------------------------------
   process an input script
------------------------------------------------------------------------- */
void Input::process_input_script()
{
    //output()->write_screen_one("hello");

    //open files for input

    {
        output().write_screen_one("\nParScale is processing your input...\n");

        //input file
        //if (inflag == 0) infile = stdin;
        //else infile = fopen(arg[inflag],"r");
        if(!infileSet_)
            infile_ = &cin;

        if (!infile_->good())
          output().write_screen_one("Cannot open input script");

        // read each line of the file

        while (!infile_->eof())
        {
            // read an entire line into memory

            char buf[MAX_CHARS_PER_LINE];
            infile_->getline(buf, MAX_CHARS_PER_LINE);

            // strip any # comment by replacing it with 0
            // do not strip # inside single/double quotes
            char *ptr = &(buf[0]);
            while (*ptr) {
              if (*ptr == '#' ) {
                *ptr = '\0';
                break;
              }
              ptr++;
            }

            // array to store memory addresses of the tokens in buf
            const char* token[MAX_TOKENS_PER_LINE] = {}; // initialize to 0

            // parse the line
            // parse the line into blank-delimited tokens
            int n = 0; // a for-loop index

            token[0] = strtok(buf, " \t"); // first token
            if (token[0]) // zero if line is blank
            {
              for (n = 1; n < MAX_TOKENS_PER_LINE; n++)
              {
                token[n] = strtok(0, " \t"); // subsequent tokens
                if (!token[n]) break; // no more tokens
              }

              //TODO error if n < 2

              // process the tokens
              std::string str(token[0], strlen(token[0]));
              pascalInterface().parse_command(str,n-1,&(token[1]));
            }


//            for (int i = 0; i < n; i++) // n = #of tokens
//              output().write_screen_one(token[i]);
          }

          output().write_screen_one("**Processed input");
          
    }
    
}

/* ----------------------------------------------------------------------
   process an input script
------------------------------------------------------------------------- */
void Input::process_JSONInput()
{
    //Open main json input file
    openJsonFile("settings", "parscale", "input", global_settings_ );
    if(global_settings_["writeJSON"].isNull())
    {
        printf("file 'settings/parscale.json' does not contain settings for input/%s. Will Abort. \n", 
                "writeJSON");
        comm().abort_one();
    }
    writeJSON_   = global_settings_["writeJSON"].toBool();    

#ifdef H5_LIB
    if(global_settings_["writeHDF5"].isNull())
    {
        printf("file 'settings/parscale.json' does not contain settings for input/%s. Will Abort. \n", 
                "writeHDF5");
        comm().abort_one();
    }
    writeHDF5_   = global_settings_["writeHDF5"].toBool();    
#else
    writeHDF5_ = false;
#endif

}

/* ----------------------------------------------------------------------
   Write JSON Containers
------------------------------------------------------------------------- */

void Input::write_containersJSON(ContainerBase &container) const
{

    if(!writeJSON_)
        return;

    const char *id = container.prop().id();
    const char *scope = container.prop().scope();
    const bool element_property = container.prop().element_property();

    // if scope equals ID, try to open JSON file here directly
    if(strcmp(id,scope) == 0)
    {

        if(container.nVec()>1)
            error().throw_error_one(FLERR,"Cannot write if container.nVec()>1.");

        //generate filenames
        char filepath[100],jsonfile[200];
        sprintf(filepath,"./%s/%.6f",
                runDirectory_,
                control().simulationState().time()
               );

        if( comm().is_proc_0() )
        {
            mkdir(filepath, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IWOTH);
            sprintf(jsonfile,"./%s/%.6f/%s.json",
                runDirectory_,
                control().simulationState().time(),
                scope
               );
        }

        comm().wait();

        if(!comm().is_proc_0())
             sprintf(jsonfile,"./%s/%.6f/%s.%d.json",
                    runDirectory_,
                    control().simulationState().time(),
                    scope,comm().me()
                    );

        QJsonObject dataObject;
        writeQJsonObject(dataObject, container);
        writeQJsonArray(jsonfile,dataObject,
                        container);

    }
}

/* ----------------------------------------------------------------------
   process JSON file
------------------------------------------------------------------------- */


void Input::openJsonFile(const char* dirName, const char* fileName, const char* objectName, QJsonObject &json ) const
{
    char jsonfile[200],errmsg[300];
    sprintf(jsonfile,"./%s/%s/%s.json",runDirectory_, dirName,fileName);

    QFile    loadFile( jsonfile );
    if(!loadFile.open(QIODevice::ReadOnly))
    {
        printf("can not open loadfile: %s. Will Abort. \n", 
                jsonfile);
        comm().abort_one();
    }

    QByteArray    saveData = loadFile.readAll();
    QJsonDocument loadDoc  = QJsonDocument::fromJson(saveData);

    if(loadDoc.isNull())
	{
        printf("QJsonDocument is invalid. Check!  \n");
        comm().abort_one();
	}

	json = loadDoc.object()[objectName].toObject();

}

// *************************************************************

void Input::writeQJsonObject(QJsonObject &json, ContainerBase &container) const
{
    json["name"]  = container.prop().scope();
    QString  str; str.setNum(control().simulationState().time(), 'g', writePrecision_);
    json["time"]  = str;
}

// * * * * * * * * * * * * * *  * * * * * * *  * * * * * * *  * * * * * * *  
void Input::writeQJsonArray(char *jsonfile,
                            QJsonObject &myParticles,
                            ContainerBase &container) const
{
    if(!container.prop().needsOutput())
        return;

    if(container.size()==0)
        return;
        
    //output().write_screen_one("Input::writeQJsonArray is writing...");
    
    //Iterate through the data for each particle
    int particleId    = 0;

    QFile saveFile( jsonfile );

    if (!saveFile.open(QIODevice::WriteOnly))
         error().throw_error_one(FLERR,"can not open savefile ", jsonfile);

    //Get access to container content
    void *ptr = container.begin_slow_dirty();   //ptr to scalar per-particle data
    QJsonObject particleDataQ;
    double doub_value;
    int    int_value;

    for(uint iPart = 0; iPart < container.size(); iPart++)
    {
        QJsonArray  dataArray;
        for(uint jVec = 0; jVec < container.lenVecUsed(); jVec++)
        {
            if(container.lenVec()>1) //these is ContainerCVODE data
            {
                if(container.isDoubleData())
                {
                    doub_value = static_cast<double***>(ptr)[iPart][0][jVec];
                    dataArray << doub_value;
                }
                else if(container.isIntData())
                {
                    int_value = static_cast<int***>(ptr)[iPart][0][jVec];
                    dataArray << int_value;
                }
                else
                        error().throw_error_one(FLERR,"Data format of container unknown.\n");
            }
            else
            {
                if(container.isDoubleData())
                {
                    doub_value = static_cast<double*>(ptr)[iPart];
                    dataArray << doub_value;                             
                }
                else if(container.isIntData())
                {
                    int_value = static_cast<int*>(ptr)[iPart];
                    dataArray << int_value;     
                }
                else
                    error().throw_error_one(FLERR,"Data format of container unknown.\n");

            }
        }

        uint globalID = particleData().returnId(iPart);
        QString  pId;  pId.setNum(globalID); 
        particleDataQ[pId] = dataArray;
    }
    myParticles["data"] = particleDataQ;

    QJsonDocument saveDoc(myParticles);
    saveFile.write( saveDoc.toJson() );
    saveFile.close();

}


// *************************************************************
void Input::readQJsonArray(const QJsonObject &myParticles, ContainerBase &container) const
{
    //Iterate through the data for each particle
    bool foundThisParticle = false;
    int particleId    = 0;
    QJsonObject::const_iterator i;
    int currentKey;

    double *  ptr     = (double*)      container.begin_slow_dirty();   //ptr to scalar per-particle data
    double ***ptr3    = (double***)    container.begin_slow_dirty();   //ptr to vectorial per-particle data

    for (i = myParticles.begin(); i != myParticles.end(); ++i)
    {
        if(particleId>=container.size())
             error().throw_error_one(FLERR,"particleId>=container.size(). You have specified too many particles in the JSON file.\n");

        currentKey    = i.key().toInt() - 1; //data indexing starts with zero
        if( currentKey>=container.size() )
        {
            printf("problem with data having the key %d: it is larger than the container.size(%d).\n", 
                   currentKey, container.size());
            error().throw_error_one(FLERR,"You have specified an incorrect key in the JSON file.\n");
        }

        QJsonArray  data = (*i).toArray();
        QJsonValue  value;
        if(data.size()<container.lenVecUsed())
        {
            printf("problem with data.size(): %d, and container.lenVecUsed(): %d.\n", 
                    data.size(), container.lenVecUsed());
            error().throw_error_one(FLERR,"j<container.lenVecUsed(). You have specified too few data in the JSON file.\n");
        }
        if(data.size()>container.lenVecUsed())
            error().throw_error_one(FLERR,"data.size()>container.lenVecUsed(). You have specified too many data in the JSON file.\n");
        for(int j = 0; j < data.size(); j++)
        {
            value = data[j];

            if(container.lenVecUsed()>1)
                 ptr3[currentKey][0][j] = value.toDouble();

            else
                 ptr[currentKey] = value.toDouble();

        }
        particleId++;
        foundThisParticle = false;
    }
    if(particleId < container.size())
        error().throw_error_one(FLERR,"particleId<container.size(). You have specified too few data in the input.");

}


// *************************************************************

void Input::fill_container_from_json(ContainerBase &container) const
{

    const char *id = container.prop().id();
    const char *scope = container.prop().scope();
    const bool element_property = container.prop().element_property();

    printf("Input::fill_container_from_json for scope: %s \n", scope);

    // if scope equals ID, try to open JSON file here directly
    if(strcmp(id,scope) == 0)
    {
        char jsonfile[200],errmsg[300];
        sprintf(jsonfile,"./%s/0/%s.json",runDirectory_, scope);

       QFile    loadFile( jsonfile );
       if(!loadFile.open(QIODevice::ReadOnly))
             error().throw_error_one(FLERR,"can not open loadfile ",jsonfile);

       QByteArray        saveData = loadFile.readAll();
       loadDoc  = QJsonDocument::fromJson(saveData);

       if(loadDoc.isNull())
             error().throw_error_one(FLERR,"QJsonDocument is invalid. Check! \n");

       myName = loadDoc.object()["name"].toString();

    }
    else
        error().throw_error_one(FLERR,"Could not read property",id,
                                "because the following JSON file was not found: ",scope);

    if(element_property)
    {
        QJsonObject myParticleData = loadDoc.object()["data"].toObject();
        if(container.isDoubleData())
        {
            readQJsonArray(myParticleData, container);
        }
        else error().throw_error_one(FLERR,"Missing implementation for",id,", ",scope);
    }
    else error().throw_error_one(FLERR,"Missing implementation for",id,", ",scope);

}

// *************************************************************
void Input::global_properties_json() const
{
	printf("\nInput::read global properties: \n");

    char global_jsonfile[200],errmsg[300];
    sprintf(global_jsonfile,"./%s/settings/parscale.json",runDirectory_);

    QFile    loadFile( global_jsonfile );
    if(!loadFile.open(QIODevice::ReadOnly))
	{
             error().throw_error_one(FLERR,"can not open loadfile ",global_jsonfile);
	}

    QByteArray        saveData = loadFile.readAll();
    loadDoc  = QJsonDocument::fromJson(saveData);

    if(loadDoc.isNull())
	{
             error().throw_error_one(FLERR,"QJsonDocument is invalid. Check! \n");
	}

    myName = loadDoc.object()["name"].toString();

	QJsonObject globalData = loadDoc.object()["global_parameters"].toObject();

	readglobal_properties(globalData);
}

// *************************************************************

void Input::readglobal_properties(const QJsonObject &json) const
{
	double global_property_1;
	global_property_1 = json["test"].toDouble();
}


//HDF5 Routines
// *************************************************************

#ifdef H5_LIB
    using namespace H5;
    #include "input_HDF5_I.h"
#else
int Input::write_containersHDF5(ContainerBase &container) const
{
    return 0;
}

#endif

