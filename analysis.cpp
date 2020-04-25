#include "includes.hpp"
#include "universe.h"
#include "analysis.h"

double Analysis::chi2teste(double *x, Universe worldTest){return -3*x[0]+6*x[0]*x[0]*x[0]*x[0]+(x[0]-1)*(x[0]-1)*(x[0]-1)-3*x[1]+6*x[1]*x[1]*x[1]*x[1]+(x[1]-1)*(x[1]-1)*(x[1]-1)-3*x[2]+6*x[2]*x[2]*x[2]*x[2]+(x[2]-1)*(x[2]-1)*(x[2]-1);};

double Analysis::f_test(double *x){return x[0]*x[0]+2*x[0]*x[0]*x[0]*x[0]+3*x[0]*x[0]*x[0];}

double Analysis::square(double x){return x*x;}

double Analysis::Chi2Bugey(double *x, Universe worldTest)
{

	
	return 0;
}



double Analysis::MinimizerMC(double (Analysis::*f) (double*, Universe), double *param, double *guess, double *inf_bounds, double *sup_bounds, int num_param, double *param_step, int steps, Universe worldTest, bool *param_run)
{



	double min=20000000;
	double localmin = 10000000;
	double init[num_param];	
	double param_min[num_param];
	double param_localmin[num_param];
	for(int i=0; i< num_param; i++){init[i]=guess[i];param_min[i]=guess[i];};
	srand(time(NULL));



	for(int i=0; i<steps; i++)
	{
		double test;


		
		for(int j=0; j<num_param; j++)
		{	
			if(param_run[j]){

				do{				
						param[j] = init[j] + param_step[j]*(1.-2.*float((rand() % 2)));
				}while( (param[j] < inf_bounds[j]) or (param[j] > sup_bounds[j]) );
			};
		};	

		test = (this->*f)(param, worldTest);


		if(test<localmin)
		{
			for(int j=0; j<num_param; j++){init[j] = param[j];param_localmin[j] = param[j];};
			localmin=test;
			if(localmin<min){
				min=localmin;
					for(int j=0; j<num_param; j++)
					{
						param_min[j] = param_localmin[j];
					};
			};		
		}
		else
		{
			int trial = (rand() % 100);
			if(trial<50)
			{
				for(int j=0; j<num_param; j++){init[j] = param[j];};
				localmin=test;
			};
		};



	};

	for(int j=0; j<num_param; j++){param[j] = param_min[j];};

	return min;

}

double Analysis::DMinimizerMC(double (Analysis::*f) (double*, Universe), double *param, double *guess, double *inf_bounds, double *sup_bounds, int num_param, double *param_step, double D, double dif, int steps, Universe worldTest, bool *param_run)
{

	double min=20000000;
	double localmin = 10000000;
	double init[num_param];	
	double param_min[num_param];
	double param_localmin[num_param];
	for(int i=0; i< num_param; i++){init[i]=guess[i];param_min[i]=guess[i];};
	srand(time(NULL));

	int count=0;

//	for(int i=0; i<steps; i++)
	do{
		double test;

		count++;
		
		for(int j=0; j<num_param; j++)
		{	
	
			if(param_run[j]){	
				do{				
						param[j] = init[j] + param_step[j]*(1.-2.*float((rand() % 2)));
				}while( (param[j] < inf_bounds[j]) or (param[j] > sup_bounds[j]) );
			};
		};
		test = (this->*f)(param, worldTest);
	

		if(abs(test-D)<localmin)
		{
			for(int j=0; j<num_param; j++){init[j] = param[j];param_localmin[j] = param[j];};
			localmin=abs(test-D);
			if(localmin<min){
				min=localmin;
					for(int j=0; j<num_param; j++)
					{
						param_min[j] = param_localmin[j];
					};
			};		
		}
		else
		{
			int trial = (rand() % 100);
			if(trial<50)
			{
				for(int j=0; j<num_param; j++){init[j] = param[j];};
				localmin=abs(test-D);
			};
		};
//		printf("Chi2: %f\tLocalmin: %f\t Min: %f\t, sin²: %f\n",test ,localmin,min, pow(10,param[5]));
		if(count>1000000){return 100;};

	}while(min>dif);
	for(int j=0; j<num_param; j++){param[j] = param_min[j];};

	return min;

}



double Analysis::MinimizerMC_2(double (Analysis::*f) (double*, Universe), double *param, double *guess, double *inf_bounds, double *sup_bounds, int num_param, double *param_step, double inf, Universe worldTest, bool *param_run)
{



	double min=20000000;
	double localmin = 10000000;
	double init[num_param];	
	double param_min[num_param];
	double param_localmin[num_param];
	for(int i=0; i< num_param; i++){init[i]=guess[i];param_min[i]=guess[i];};
	srand(time(NULL));

	double test;

	do
	{



		
		for(int j=0; j<num_param; j++)
		{	
			if(param_run[j]){

				do{				
						param[j] = init[j] + param_step[j]*(1.-2.*float((rand() % 2)));
				}while( (param[j] < inf_bounds[j]) or (param[j] > sup_bounds[j]) );
			};
		};	

		test = (this->*f)(param, worldTest);


		if(test<min)
		{
			for(int j=0; j<num_param; j++){init[j] = param[j];param_localmin[j] = param[j];};
			if(test<min){
				min=test;
					for(int j=0; j<num_param; j++)
					{
						param_min[j] = param_localmin[j];
					};
			};		
		};

	}while(abs(test-min)<0.01);

	for(int j=0; j<num_param; j++){param[j] = param_min[j];};

	return min;

}




Analysis::Analysis(){}

Analysis::~Analysis(){}


