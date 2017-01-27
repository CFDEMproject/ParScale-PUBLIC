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


#include "coupling_model.h"
#include "string.h"
#include "comm.h"
#include "input.h"
#include "output.h"
#include "error.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

CouplingModel::CouplingModel(ParScale *ptr,const char *name) :
    ParScaleBase(ptr),
    name_(0),
    scope_(0),
    verbose_(false),
    pullEvery_(1),
    pushEvery_(1),
    pullNext_(1),
    hasPulled_(false),
    pushNext_(1),
    hasPushed_(false),
    exchangeEventsLocalId_(NULL),
    exchangeEventsReceivingProcess_(NULL)
{

    int n = strlen(name) + 1;
    name_ = new char[n];
    strcpy(name_,name);

    n = strlen(name) + 9 + 1;
    scope_ = new char[n];
    strcpy(scope_,"coupling_");
    strcat(scope_,name);


    char msgstr[500];
    sprintf(msgstr,"CouplingModel with name %s initialized. \n", name_);
    output().write_screen_all(msgstr);

//    if(comm().is_parallel())
//        error().throw_error_all(FLERR,"TODO: need to communicate settings in parallel in Input:: function");
}

/* --------------------------------------------------------------------- */
CouplingModel::~CouplingModel()
{
    delete [] name_;
    delete [] scope_;
}

/* --------------------------------------------------------------------- */
void CouplingModel::read()
{
    //common properties
    input().openJsonFile("settings", scope(), "properties", properties_ );
    verbose_   = properties_["verbose"].toBool();
    pullEvery_ = properties_["pullEvery"].toInt();

    if ( strcmp(name_, "liggghts") == 0)    //has to be set to 1 due to the conductive flux management of liggghts
        pullEvery_ = 1;         //TODO in future: sum flux from liggghts locally in PartScale to avoid a coupling at every time step

    pushEvery_ = properties_["pushEvery"].toInt();
    output().write_screen_all("read in for coupling model done \n");

}

/* --------------------------------------------------------------------- */
void CouplingModel::setPushPull() const
{
    if(hasPulled_)
        pullNext_ += pullEvery_;

    if(hasPushed_)
        pushNext_ += pushEvery_;

    hasPulled_ = false;
    hasPushed_ = false;
}
