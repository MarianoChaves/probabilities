#pragma once

class Detector{

	private:
		vector <double> statError;
		vector <double> sysError;
		vector <double> error;

		vector<double> data;
		vector<double> energyBins;
		int Nbins;

	
	public:


		void setError(vector <double> error){this->error=error;};
		vector <double> getError(){return this->error;};

		void setSysError(vector <double> sysError){this->sysError=sysError;};
		vector <double> getSysError(){return this->sysError;};

		void setStatError(vector <double> statError){this->statError=statError;};
		vector <double> getStatError(){return this->statError;};

		void setData(vector <double> data){this->data=data;};
		vector <double> getData(){return data;};

		void setEnergyBins(vector <double> energyBins){this->energyBins=energyBins;};
		vector <double> getEnergyBins(){return this->energyBins;};

		void setNbins(int Nbins){this->Nbins=Nbins;};
		int getNbins(){return this->Nbins;};

		Detector(char *filename);
		Detector();
		~Detector();

};
