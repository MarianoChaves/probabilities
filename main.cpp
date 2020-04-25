#include "includes.hpp"
#include "universe.h"
#include "detector.h"
#include "analysis.h"



/* ############################################################## */
/* ############################################################## */




int main(int argc, char *argv[])
{ 	


/***************************************************************************
 *                        I N I T I A L I Z A T I O N                      *
 ***************************************************************************/

ofstream outfile1("plots/prob.dat");


Universe world_SO;
Universe world_NSI;

double theta12=asin(sqrt(0.31));
double theta13=asin(sqrt(0.023));
double theta23=asin(sqrt(0.5));
double deltaCP=0;
double dm21=7.5e-5;
double dm31=2.4e-3;

world_SO.setDm21(dm21);
world_SO.setDm31(dm31);
world_SO.setTh12(theta12);
world_SO.setTh13(theta13);
world_SO.setTh23(theta23);
world_SO.setDcp(deltaCP); 

double density[1] = {2.8};
double length[1] = {1300};
world_SO.setDensity(density);
world_SO.setLength(length);

world_NSI=world_SO;

double nsi_SO[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
world_SO.setNewPhysics(nsi_SO);


double ee=0.0;
double eu=0.1;
double et=0.5;
double tt=0.0;
double ut=0.0;
double nsi[5] = {ee, eu, et, tt, ut};
world_NSI.setNewPhysics(nsi);




double Pnum,Pmy, PSO, PmySO;
int N = 500;

double E_i= log10(0.04);
double E_f= log10(10.0);
double dE = (E_f-E_i)/N;

for(int i; i<N; i++)
{	
	double En=pow(10.,E_i+dE*i);
	Pnum = world_NSI.getProbability(En,1,0,1);
	Pmy = world_NSI.getMyProbability(En,1);

	PSO = world_SO.getProbability(En,1,0,1);
	PmySO = world_SO.getMyProbability(En,1);
	outfile1<<En<<"\t"<<Pnum<<"\t"<<Pmy<<"\t"<<PSO<<"\t"<<PmySO<<endl;
};

outfile1.close();

}
// main program ends!

/* ############################################################## */
/* ############################################################## */


