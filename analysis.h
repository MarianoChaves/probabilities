#include "detector.h"

class Analysis{

	private:

		Detector *exp;

		
	
	public:

		double MinimizerMC(double (Analysis::*f)(double*, Universe), double *param, double *guess, double *inf_bounds, double *sup_bounds, int num_param, double *param_step, int steps, Universe worldTest, bool *param_run);

		double MinimizerMC_2(double (Analysis::*f)(double*, Universe), double *param, double *guess, double *inf_bounds, double *sup_bounds, int num_param, double *param_step, double inf, Universe worldTest, bool *param_run);

		double DMinimizerMC(double (Analysis::*f)(double*, Universe), double *param, double *guess, double *inf_bounds, double *sup_bounds, int num_param, double *param_step, double D, double CL, int steps, Universe worldTest, bool *param_run);


		double chi2teste(double *param, Universe worldTest);

		double Chi2Bugey(double *param, Universe worldTest);

		double DChi2Bugey(double *param, double D , Universe worldTest);

		void setExp(Detector *exp){this->exp=exp;};
		Detector* getExp(){return this->exp;};

	
		double square(double x);
		double f_test(double *x);

		Analysis();
		~Analysis();

};
