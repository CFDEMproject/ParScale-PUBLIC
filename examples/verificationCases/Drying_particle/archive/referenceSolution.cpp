#include<stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace std;
      
    FILE *fp;  
      
// *** USER INPUT *************************************************
// Physical Parameters
    const double t_end  =   0.3;     
    const double R      =   5e-3;    //sphere radius
    const double alpha  =   500.0;   //heat transfer coefficient
    const double lambda =   5.0;
    const double rhoP   =   1000;
    const double cpP    =   300.0;
    const double temp_sphere        = 800.0; 
    const double temp_environment   = 300.0; 

// Numerical Parameters
    const int    size_2 = 25 ;       //radial grid points
    const int    id2Plot= 14;        //radial index to plot

// *** USER INPUT - END ********************************************

    double size_1;
    double r;
    double temp[2][size_2+1];
    double dt;
    double Fo;
    
    //calculation of constants
    const double a = lambda/rhoP/cpP;   //thermal diffusivity, a = lambda/rho/cp
    const double dr=(R/size_2);
	const double Bi=(alpha*dr)/lambda;
	const double dt_1 = ((dr*dr)/(6.*a));
    const double dt_3 = (((1.*R*dr)/(a*2.*Bi))/(((R*dr)/a)*((1.*a/(Bi*dr*dr))-(a/(dr*dr)))+1))-1e-06;  
    //double dt_2 =((dr*dr)/(2.*a))-1e-06 - immer groesser als dt_1 & dt_3;

    //variables
    int j; //radius
    int t; //time

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  	
int write(double tdt)
{
    fp = fopen("referenceSolution.dat","a");
    fprintf(fp,"r\t T(r) for time =  \n", tdt);
    for (j=0; j<size_2+1; j++)
    	fprintf(fp,"%lf \t %lf\n",(j*dr)/R,temp[1][j]);
    fclose(fp);
    
    fp = fopen("referenceSolution.center.dat","a");
   	fprintf(fp,"%g \t %g \n", tdt,temp[1][0]);
    fclose(fp);

    fp = fopen("referenceSolution.inBetween.dat","a");
   	fprintf(fp,"%g \t %g \n", tdt,temp[1][id2Plot]);
    fclose(fp);

    fp = fopen("referenceSolution.surface.dat","a");
   	fprintf(fp,"%g \t %g \n", tdt,temp[1][size_2]);
    fclose(fp);
    
    return(0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
int main(void)
{
    if (dt_1<=dt_3)                                            
       dt=dt_1;
    else
      dt=dt_3;

      //calculation Fouier number and number of time steps
    Fo=((a*dt)/(dr*dr));
    size_1=(t_end/dt);

    for (j=0; j<size_2+1; j++)                                      
    {
		temp[0][j]= (temp_sphere);
	    temp[1][j]= (temp_sphere);
	}
      
    for (t=0; t<size_1+1; t++)                                      
    {
        //explicite solution at r=0 with bc 2nd order
        temp[1][0]=(temp[0][0]*(-6.*Fo+1)+6.*Fo*temp[0][1]);         

        for (j=1; j<size_2; j++)                                    
        {
    		r=j*dr;
    		
    		//explicite solution for all inner knots
            temp[1][j]=( temp[0][j]*(1.-2.*Fo)
                         +Fo*(temp[0][j+1]+temp[0][j-1])
                         +a*dt*((temp[0][j+1]-temp[0][j-1])/(r*dr))
                       );                
        }              
        //explicite solution at r=R (outer knot) with bc 3rd order
		//temp_environment = temp_environment_begin + ((temp_environment_end-temp_environment_begin)/size_1)*t;
		//printf("actual temp water = %g ", temp_environment);
        temp[1][size_2]=  temp[0][size_2]*
                          (1.-2.*Fo-2.*Fo*Bi-2.*Bi*((a*dt)/(R*dr)))
                          +Fo*(2.*Bi*temp_environment+2.*temp[0][size_2-1])
                          +((a*dt)/(R*dr))*2.*Bi*temp_environment;  

        write(t*dt);
		
	    for (j=0; j<size_2+1; j++)                                  
 		    temp[0][j]=temp[1][j];
    }
    
	return(0);
}
