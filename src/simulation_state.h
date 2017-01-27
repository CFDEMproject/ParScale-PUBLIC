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

#ifndef PASC_SIMULATION_STATE_H
#define PASC_SIMULATION_STATE_H

#include "stdio.h"
#include "pascal_base.h"

namespace PASCAL_NS
{

class SimulationState : public ParScaleBase
{
    public:

      SimulationState(ParScale *ptr);

      void setDeltaT(double _dt)
      { deltaT_ = _dt; }


      void setOutputTimeStep(double _dt)
      {
         deltaT_writeContainers_ = _dt;
         time_writeContainers_   = time_ + _dt ;
      }


      void progress_time(double _dt)
      { time_ += _dt; time_step_++; }

      inline double time()
      { return time_; }

      inline double deltaT()
      { return deltaT_; }

      inline int timeStep()
      { return time_step_; }

      inline bool writeContainers()
      {
        if(time_ >= time_writeContainers_)
        {
            time_writeContainers_ += deltaT_writeContainers_;
            return true;
        }

        return false;
      }


    private:

      double time_;        //the current time
      double deltaT_;      //the time step
      int    time_step_;   //the current integration step

      //I/O information
      double deltaT_writeContainers_;     //the time intervall for writing containers
      double time_writeContainers_;       //the next time for writing containers


};

} //end PASCAL_NS

#endif
