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

#ifndef PASC_INPUT_H
#define PASC_INPUT_H


#include "stdio.h"
#include "error.h"
#include "input_base.h"
#include "pascal_base.h"
#include "container.h"
#include "pascal_base_interface.h"
#include <fstream>

#ifdef H5_LIB
    #include "hdf5.h"
    #include "H5Cpp.h"
#endif

#define __PIC__
#include <QtCore/QJsonArray>
#include <QtCore/QJsonObject>
#include <QtCore/QJsonDocument>
#include <QtCore/QString>
#include <QtCore/QFile>

using std::ifstream;

namespace PASCAL_NS
{

class Input : public ParScaleBase, public ParScaleBaseInterface, public InputBase
{
    public:

      Input(ParScale *ptr);
      ~Input();

      void set_input_script(const char *);  //sets the name of the input script

      void set_runDirectory(const char *);  //sets the name of the directory to run in

      char* runDirectory() {return runDirectory_;};

      void process_input_script();

      void process_input(const char *);  //processes a single input command

      void process_JSONInput();

      // read a specified input type from JSON format and store it
      // in a container. storage in a container is preferred so can
      // comunicate, restart more easily
      void fill_container_from_json(class ContainerBase &container) const;
      
      void write_containersJSON(class ContainerBase &container) const;

      void openJsonFile(const char*, const char*, const char*, QJsonObject &json ) const;

      void readQJsonArray(const QJsonObject &json, ContainerBase &container) const;

      void writeQJsonObject(QJsonObject &json, ContainerBase &container) const;

      void writeQJsonArray(char *jsonfile, QJsonObject &json, 
                           ContainerBase &container) const;

	  void global_properties_json() const;
	 
      void readglobal_properties(const QJsonObject &json) const;

      //under namespace H5
      int write_containersHDF5(class ContainerBase &container) const;


    private:

      char   *runDirectory_;            // (relative) directory path to run in

      bool    infileSet_;

      ifstream myFile;                  // infile - file stream

      istream *infile_;                 // input stream

      mutable QJsonDocument  loadDoc;

      QJsonObject            global_settings_;
      
      bool                   writeJSON_;
      
      bool                   writeHDF5_;

      mutable QString        myName;
      
      mutable int            writePrecision_;
};

  // *************************************
  #include "input_I.h"
  // *************************************

} //end PASCAL_NS

#endif
