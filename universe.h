#pragma once

class Universe{
	
	private:
		double dm31, dm21;	
		double th13, th12, th23;
		double dcp;
		complex <double> U[3][3];
		double P[3][3];

		double ee, eu, et, tt, ut;

		double *density;	
		double *length;
		int psteps;

	public:
		void setDm31(double dm31){this->dm31=dm31;};
		void setDm21(double dm21){this->dm21=dm21;};
		void setTh12(double th12){this->th12=th12;};
		void setTh13(double th13){this->th13=th13;};
		void setTh23(double th23){this->th23=th23;};
		void setDcp(double dcp){this->dcp=dcp;};

		void setNewPhysics(double *newp){this->ee=newp[0];this->eu=newp[1];this->et=newp[2];this->tt=newp[3];this->ut=newp[4];}
		
		void setDensity(double *density){this->density=density;};
		void setLength(double *length){this->length=length;this->psteps=sizeof(length);};

		double getProbability(double E, int cp_sign, int nu_i, int nu_f);

		double getMyProbability(double E,int cp_sign);	

		double square(double x);
		double min(double x, double y);	
		
		Universe();
		~Universe();

}; 
