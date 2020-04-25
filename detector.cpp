#include "includes.hpp"
#include "universe.h"
#include "detector.h"



Detector::Detector(){}

Detector::Detector(char *filename)
{
	FILE *datafile;
	datafile = fopen(filename,"r");

	fscanf(datafile,"number of bins: %d\n",&Nbins);
	float aux;

	fscanf(datafile,"data: {");
	for(int i = 0; i< Nbins; i++)
	{
		fscanf(datafile,"%f, ", &aux);
		data.push_back(double(aux));

	};

	fscanf(datafile,"}\n\ncenter of energy bins: {");
	for(int i = 0; i< Nbins; i++)
	{
		fscanf(datafile,"%f, ", &aux);
		energyBins.push_back(double(aux));
	};


	fscanf(datafile,"}\n\nsystematic errors: {");
	for(int i = 0; i< Nbins; i++)
	{
		fscanf(datafile,"%f, ", &aux);
		sysError.push_back(double(aux));

	};

	fscanf(datafile,"}\n\nstatistical errors: {");
	for(int i = 0; i< Nbins; i++)
	{
		fscanf(datafile,"%f, ", &aux);
		statError.push_back(double(aux));
	};

	for(int i = 0; i< Nbins; i++)
	{		
		error.push_back(sqrt(statError[i]*statError[i] + sysError[i]*sysError[i]));

	};

	fclose(datafile);
}

Detector::~Detector(){}

