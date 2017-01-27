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



#include "model_base.h"
#include "string.h"
#include "comm.h"
#include "output.h"
#include "input.h"
#include "error.h"
#include <fstream>

using namespace PASCAL_NS;
using std::ifstream;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

ModelBase::ModelBase(ParScale *ptr,const char *name) : ParScaleBaseAccessible(ptr),
    name_(0),
    scope_(0),
    alwaysTrueVariable_(true)
{
    int n = strlen(name) + 1;
    name_ = new char[n];
    strcpy(name_,name);

    n = strlen(name) + 6 + 1;
    scope_ = new char[n];
    strcpy(scope_,"model_");
    strcat(scope_,name);

    char msgstr[500];
    sprintf(msgstr,"Model with name %s initialized. \n", name_);
    output().write_screen_all(msgstr);

//    if(comm().is_parallel())
//        error().throw_error_all(FLERR,"TODO: need to communicate settings in parallel in Input:: function");
}

// ----------------------------------------------------------------------
void ModelBase::read_model_json_file(const char *model_name, double *ptr, vector<double> &parameters)
{
    //TODO: this is ugly, do similar to  input().openJsonFile

    char jsonfile[200];
    sprintf(jsonfile,"./%s/settings/model_%s.json", input().runDirectory(), model_name);

    QFile    loadFile( jsonfile );
    if(!loadFile.open(QIODevice::ReadOnly))
             error().throw_error_one(FLERR,"can not open loadfile ",jsonfile);

    QByteArray  saveData = loadFile.readAll();
                loadDoc  = QJsonDocument::fromJson(saveData);

    if(loadDoc.isNull())
    {
        error().throw_error_one(FLERR,"QJsonDocument is invalid. Check! \n");
    }
    /*else
    {
        error().throw_error_one(FLERR,"Could not read constants for model ",name_model,
                                "because the following JSON file was not found: ",jsonfile);
    }*/
    readQJsonObject(loadDoc.object());
    QJsonObject constant_model_properties = loadDoc.object()["properties"].toObject();
    readQJsonConstant(constant_model_properties, model_name,ptr);
    readQJsonParameterVector(constant_model_properties, model_name, parameters);

}

// ----------------------------------------------------------------------
QJsonObject ModelBase::readQJsonObject(const char* _file_name, const char *_object_name)
{
    char jsonfile[200];
    sprintf(jsonfile,"./%s/settings/%s.json", input().runDirectory(), _file_name);

    QFile    loadFile( jsonfile );
    if(!loadFile.open(QIODevice::ReadOnly))
             error().throw_error_one(FLERR,"Please supply this file:",jsonfile);

   QByteArray  saveData = loadFile.readAll();
               loadDoc  = QJsonDocument::fromJson(saveData);

    if(loadDoc.isNull())
        error().throw_error_one(FLERR,"QJsonDocument is invalid. Check this file:", jsonfile);

    readQJsonObject(loadDoc.object());

    QJsonObject myOb = loadDoc.object()[_object_name].toObject();

    if(myOb.isEmpty())
        error().throw_error_one(FLERR,"QJsonObject is empty. Check to be sure to have this object in the file:",
                                _object_name,
                                jsonfile);
    return myOb;
}



// ----------------------------------------------------------------------
void ModelBase::read_verbose_json_file(const char *property_name, bool *ptr)
{

    QJsonObject property_value_ = readQJsonObject("verbose", "verbose");
    if(property_value_[property_name].isNull())
      error().throw_error_one(FLERR,"ERROR: property_name not found in file. \n",
                              property_name,
                              "settings/verbose.json");

    QJvalue_ = property_value_[property_name].toBool();
    *ptr=QJvalue_.toBool();
}

// ----------------------------------------------------------------------
void ModelBase::read_chemistry_single_react_json_file(const char *property_name, double *ptr, bool strict)
{

    QJsonObject property_value_ = readQJsonObject("chemistry_single_reaction", "reaction");
    if(property_value_[property_name].isNull() && strict)
      error().throw_error_one(FLERR,"ERROR: property_name not found in file. \n",
                              property_name,
                              "settings/chemistry_single_reaction.json");

    else if(property_value_[property_name].isNull()) //return if not specified but non-strict
        return;

    QJvalue_ = property_value_[property_name].toDouble();
    *ptr=QJvalue_.toDouble();
}

// ----------------------------------------------------------------------
void ModelBase::read_chemistry_single_react_json_file(const char *property_name, bool *ptr, bool strict)
{
    QJsonObject property_value_ = readQJsonObject("chemistry_single_reaction", "reaction");
    if(property_value_[property_name].isNull() && strict)
      error().throw_error_one(FLERR,"ERROR: property_name not found in file. \n",
                              property_name,
                              "settings/chemistry_single_reaction.json");

    else if(property_value_[property_name].isNull())  //return if not specified but non-strict
        return;

    QJvalue_                   = property_value_[property_name].toBool();
    QJsonValue  QJvalueDouble_ = property_value_[property_name].toDouble(); //if set to a value, check if positive
    if(QJvalueDouble_.toDouble()>0)
        *ptr = &(alwaysTrueVariable_);
    else
        *ptr=QJvalue_.toBool();
}


// ----------------------------------------------------------------------
void ModelBase::read_species_properties(const char *species_name, const char *property_name, double *ptr)
{

    QJsonObject property_value_ = readQJsonObject("species_properties", species_name);

    if(property_value_[property_name].isNull())
      error().throw_error_one(FLERR,"ERROR for a species: property_name not found in file. \n",
                              species_name,
                              property_name,
                              "settings/species_properties.json");

    QJvalue_ = property_value_[property_name].toDouble();
    *ptr=QJvalue_.toDouble();
}

// ----------------------------------------------------------------------
void ModelBase::readQJsonConstant(const QJsonObject &json, const char *name_model,double *ptr)
{
    model_constant_ = json["model"].toString();
    QString qt_const_ = "const";

    if (strcmp(qPrintable(model_constant_),qPrintable(qt_const_)) == 0)
    {
        QJvalue_ = json["const_value"].toDouble();
        *ptr=QJvalue_.toDouble();
        printf("Assuming constant model for %s with default value of: %g \n",name_model, *ptr);
    }

    else if (json.size()==0)
    {
        printf("Assuming %s is global property file?\n",name_model);
    }

    else if (strcmp(qPrintable(model_constant_),qPrintable(qt_const_)) != 0 && (json.size()==2))
    {
        printf("Non-constant properties not supported yet. Please check for model:const in your json-file %s \n",name_model);
        error().throw_error_one(FLERR,"ERROR in json-file \n",name_model);
    }

    //printf("json.size = %i \n ",json.size());
}

// ----------------------------------------------------------------------
void ModelBase::readQJsonParameterVector(const QJsonObject &json, const char *name_model, vector<double> &parameters)
{
    QJsonArray  myParams = json["parameters"].toArray();

    QJsonArray::const_iterator i;
    for (i = myParams.begin(); i != myParams.end(); ++i)
    {
        parameters.push_back((*i).toDouble());
    }
}


// ----------------------------------------------------------------------
void ModelBase::readQJsonObject(const QJsonObject &json)
{
    myName = json["name"].toString();
    //printf("name: %s \n", qPrintable(myName));
}

// ----------------------------------------------------------------------
ModelBase::~ModelBase()
{
    delete [] name_;
    delete [] scope_;
}
