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
#ifdef COUPLING_MODEL_CLASS

CouplingModelStyle(json, CouplingModelJSON)

#else


#ifndef PASC_COUPLING_MODEL_JSON_H
#define PASC_COUPLING_MODEL_JSON_H

#include "coupling_model.h"
#include "qjson_includes.h"

namespace PASCAL_NS
{

class CouplingModelJSON : public CouplingModel
{
    public:

      CouplingModelJSON(ParScale *ptr, const char *_name);

      void    read();

      virtual bool external_code_in_control() const
      { return false; }

    private:

      // pull
      bool fill_container_from_coupling(class ContainerBase &container) const;

      // push
      bool dump_container_to_coupling(class ContainerBase &container) const;

      // force a fill operation for all containers,  MUST DO for JSON-type coupling!
      bool forceFill() const {return true;};

      void insertSingleValue(double currValue, class ContainerBase &container, int, double*)  const;

      void insertSingleValueInArray(double currValue, class ContainerBase &container, int, double***)  const;

      //settings for BCs
      mutable int applyAt_;
      QJsonObject setBCs_;
      vector<QString>      settingBCName_;
      vector<QString>      settingBCType_;
      vector<int>          settingBCApplyAt_;
      vector<QString>      settingBCEqn_;
      vector<QJsonArray>   settingBCTime_;    //holds vector with times at which settings are applied
      vector<QJsonArray>   settingBCValue_;   //holds vector with values to apply
};

} //end PASCAL_NS

#endif

#endif
