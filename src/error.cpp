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


#include "error.h"
#include "input.h"
#include "output.h"
#include "comm.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

Error::Error(ParScale *ptr) : ParScaleBase(ptr)
{

}

/* ----------------------------------------------------------------------
   called by one proc
------------------------------------------------------------------------- */

void Error::throw_error_one(const char *file, int line,const char *msg1,
                            const char *msg2,const char *msg3,const char *msg4) const
{
    char errstr[500];
    sprintf(errstr,"\n\nERROR on process %d (in file %s, line %d): ",comm().me(),file,line);

    // just one proc is calling, so call all() write
    output().write_screen_all(errstr);
    output().write_screen_all(msg1);
    if(msg2) output().write_screen_all(msg2);
    if(msg3) output().write_screen_all(msg3);
    if(msg4) output().write_screen_all(msg4);
    output().write_log_all(errstr);
    output().write_log_all(msg1);
    if(msg2) output().write_log_all(msg2);
    if(msg3) output().write_log_all(msg3);
    if(msg4) output().write_log_all(msg4);

    int a = comm().nprocs();
    comm().abort_one();
}

/* ----------------------------------------------------------------------
   called by all proc
------------------------------------------------------------------------- */

void Error::throw_error_all(const char *file, int line,const char *msg1,
                            const char *msg2,const char *msg3,const char *msg4) const
{
    char errstr[500];
    sprintf(errstr,"\n\nERROR on process %d (in file %s, line %d): ",comm().me(),file,line);

    // all procs are calling, just 0 should write
    output().write_screen_one(errstr);
    output().write_screen_one(msg1);
    if(msg2) output().write_screen_all(msg2);
    if(msg3) output().write_screen_all(msg3);
    if(msg4) output().write_screen_all(msg4);
    output().write_log_one(errstr);
    output().write_log_one(msg1);
    if(msg2) output().write_log_all(msg2);
    if(msg3) output().write_log_all(msg3);
    if(msg4) output().write_log_all(msg4);

    comm().finalize_all();
}
