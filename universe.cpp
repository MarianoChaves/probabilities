#include "includes.hpp"
#include "universe.h"
#include "aux/zheevq3.h"

#define EXP_15 0
#define EXP_45 1
#define EXP_95 2


#define SQR(x)      ((x)*(x))                        // x^2 
#define SQR_ABS(x)  (SQR(real(x)) + SQR(imag(x)))  // |x|^2

Universe::Universe(){}

Universe::~Universe(){}

double Universe::getMyProbability(double E, int cp_sign)
{

	complex<double> I{0.0, 1.0};
	double V = cp_sign*0.475*0.000076*density[0];
	double A = 2*E*V;
	double L = length[0]*1.27*4;
	double s13=sin(th13);
	double c23=cos(th23);
	double s23=sin(th23);
	double c12=cos(th12);
	double s12=sin(th12);

	double rD= dm21/dm31;
	double rA= A/dm31;
	double D31 = dm31/(2*E);

	double etT=c23*et;
	double euT=-s23*et;
	double ttT=c23*c23*tt+2*c23*s23*real(ut);
	double utT=-s23*c23*tt+(c23*c23-s23*s23)*ut;
	//double utT=-s23*c23*tt+(c23*c23-s23*s23)*real(ut)+I*imag(ut);

	double Sigma = abs(s13*exp(-I*dcp)+rA*etT);
	double Omega = abs(rD*c12*s12+rA*euT);
	double Lambda = 1./rA+ttT;
	double Gamma = 1+ee;
	double eta = Lambda-Gamma;

	double Sigma_cp = arg(s13*exp(-I*dcp)+rA*etT);
	double Omega_cp = arg(rD*c12*s12+rA*euT);
	double utT_cp = arg(utT);
	
	double P, P1, P32, P2, P2_1, P2_2, P2_3;

	P1 = 4*square(Sigma*s23/(rA*eta)*sin(D31*L*rA*eta/2));

	P32 = 8*c23*s23*Sigma*Omega*sin(D31*L*rA*Gamma/2)*sin(D31*L*rA*eta/2)*cos(D31*L*rA*Lambda/2-Sigma_cp+Omega_cp)/(rA*rA*Gamma*eta); 

	P2_1 = 4*square(c23*Omega*sin(D31*L*rA*Gamma/2)/(rA*Gamma))+2*square(Sigma*s23)*(2*square(Sigma)/pow(rA*eta,3)-(rD*s12*s12+2*s13*s13)/square(rA*eta))*(D31*L)*(sin(rA*eta*D31*L));

	P2_2 = -4*s23*s23*(4*pow(Sigma/(rA*eta),4)-2*square(Sigma)*(rD*s12*s12+2*s13*s13)/pow(rA*eta,3)+Sigma*s13*(2*rD*s12*s12+s13*s13)*cos(dcp+Sigma_cp)/square(rA*eta))*square(sin(D31*L*rA*eta/2));

	P2_3 = 4*c23*utT*square(Sigma)*s23*sin(D31*L*rA*eta/2)*(sin(utT_cp-D31*L*rA*(Lambda+Gamma)/2)/(rA*rA*eta*Gamma*Lambda)-sin(utT_cp-D31*L*rA*eta/2)/(rA*rA*eta*eta*Gamma)+sin(utT_cp+D31*L*rA*eta/2)/(rA*rA*eta*eta*Lambda));

	P2=P2_1+P2_2+P2_3;
	P=P1+P32+P2;
	return P;



}

double Universe::getProbability(double E, int cp_sign, int nu_i, int nu_f)
{

	complex <double> A[3][3]; 
	complex <double> Q[3][3];
	double P[3][3];
	double w[3];
	double Ve,Vee,Veu,Vet,Vtt,Vut;
	double sin12,sin13,sin23,cos12,cos13,cos23;

	complex<double> I{0.0, 1.0};

	  int z, L, ii, jj;
	  /* Set all probabilities to zero initially */
	  for (ii=0; ii < 3; ii++)
	  {
	    for (jj=0; jj < 3; jj++)
	  {
	      P[ii][jj] = 0.0;

	  }
	  }





	complex <double> Pee[psteps];
	complex <double> Pem[psteps];
	complex <double> Pet[psteps];

	complex <double> Pme[psteps];
	complex <double> Pmm[psteps];
	complex <double> Pmt[psteps];

	complex <double> Pte[psteps];
	complex <double> Ptm[psteps];
	complex <double> Ptt[psteps];


	for (z=0; z < psteps; z++)
		{ 
		L += length[z];

		sin12=sin(th12);
		sin13=sin(th13);
		sin23=sin(th23);
		cos12=cos(th12);
		cos13=cos(th13);
		cos23=cos(th23);


		complex <double> Ue1 = cos12*cos13;
		complex <double> Ue2 = sin12*cos13;
		complex <double> Ue3 = sin13*exp(-I*double(cp_sign)*dcp);

		complex <double> Um1 = -sin12*cos23-cos12*sin23*sin13*exp(I*double(cp_sign)*dcp);
		complex <double> Um2 = cos12*cos23-sin12*sin23*sin13*exp(I*double(cp_sign)*dcp);
		complex <double> Um3 = sin23*cos13;

		complex <double> Ut1 = sin12*sin23-cos12*cos23*sin13*exp(I*double(cp_sign)*dcp);
		complex <double> Ut2 = -cos12*sin23-sin12*cos23*sin13*exp(I*double(cp_sign)*dcp);
		complex <double> Ut3 = cos23*cos13;


		Ve=cp_sign*0.475*0.000076*density[z];
		Vee=cp_sign*0.475*ee*0.000076*density[z];
		Veu=cp_sign*0.475*eu*0.000076*density[z];
		Vet=cp_sign*0.475*et*0.000076*density[z];
		Vtt=cp_sign*0.475*tt*0.000076*density[z];
		Vut=cp_sign*0.475*ut*0.000076*density[z];


		A[0][0] = (0.5/E)*(Ue2*conj(Ue2)*dm21+Ue3*conj(Ue3)*dm31) + Ve + Vee ;
		A[0][1] = (0.5/E)*(Ue2*conj(Um2)*dm21+Ue3*conj(Um3)*dm31) + Veu;
		A[0][2] = (0.5/E)*(Ue2*conj(Ut2)*dm21+Ue3*conj(Ut3)*dm31) + Vet;

		A[1][0] = (0.5/E)*(Um2*conj(Ue2)*dm21+Um3*conj(Ue3)*dm31) + Veu ;
		A[1][1] = (0.5/E)*(Um2*conj(Um2)*dm21+Um3*conj(Um3)*dm31) ;
		A[1][2] = (0.5/E)*(Um2*conj(Ut2)*dm21+Um3*conj(Ut3)*dm31) + Vut;

		A[2][0] = (0.5/E)*(Ut2*conj(Ue2)*dm21+Ut3*conj(Ue3)*dm31) + Vet;
		A[2][1] = (0.5/E)*(Ut2*conj(Um2)*dm21+Ut3*conj(Um3)*dm31) + Vut;
		A[2][2] = (0.5/E)*(Ut2*conj(Ut2)*dm21+Ut3*conj(Ut3)*dm31) + Vtt;


		zheevq3(A, Q, w);


		double L1 = w[0];
		double L2 = w[1];
		double L3 = w[2];

		complex <double> Ue1f = Q[0][0];
		complex <double> Ue2f = Q[0][1];
		complex <double> Ue3f = Q[0][2];

		complex <double> Um1f = Q[1][0];
		complex <double> Um2f = Q[1][1];
		complex <double> Um3f = Q[1][2];

		complex <double> Ut1f = Q[2][0];
		complex <double> Ut2f = Q[2][1];
		complex <double> Ut3f = Q[2][2];



		Pee[z] = conj(Ue1f)*Ue1f*exp(-I*4.*1.27*L1*length[z])+conj(Ue2f)*Ue2f*exp(-I*4.*1.27*L2*length[z])+conj(Ue3f)*Ue3f*exp(-I*4.*1.27*L3*length[z]);
		Pem[z] = conj(Ue1f)*Um1f*exp(-I*4.*1.27*L1*length[z])+conj(Ue2f)*Um2f*exp(-I*4.*1.27*L2*length[z])+conj(Ue3f)*Um3f*exp(-I*4.*1.27*L3*length[z]);
		Pet[z] = conj(Ue1f)*Ut1f*exp(-I*4.*1.27*L1*length[z])+conj(Ue2f)*Ut2f*exp(-I*4.*1.27*L2*length[z])+conj(Ue3f)*Ut3f*exp(-I*4.*1.27*L3*length[z]);

		Pme[z] = conj(Um1f)*Ue1f*exp(-I*4.*1.27*L1*length[z])+conj(Um2f)*Ue2f*exp(-I*4.*1.27*L2*length[z])+conj(Um3f)*Ue3f*exp(-I*4.*1.27*L3*length[z]);
		Pmm[z] = conj(Um1f)*Um1f*exp(-I*4.*1.27*L1*length[z])+conj(Um2f)*Um2f*exp(-I*4.*1.27*L2*length[z])+conj(Um3f)*Um3f*exp(-I*4.*1.27*L3*length[z]);
		Pmt[z] = conj(Um1f)*Ut1f*exp(-I*4.*1.27*L1*length[z])+conj(Um2f)*Ut2f*exp(-I*4.*1.27*L2*length[z])+conj(Um3f)*Ut3f*exp(-I*4.*1.27*L3*length[z]);

		Pte[z] = conj(Ut1f)*Ue1f*exp(-I*4.*1.27*L1*length[z])+conj(Ut2f)*Ue2f*exp(-I*4.*1.27*L2*length[z])+conj(Ut3f)*Ue3f*exp(-I*4.*1.27*L3*length[z]);
		Ptm[z] = conj(Ut1f)*Um1f*exp(-I*4.*1.27*L1*length[z])+conj(Ut2f)*Um2f*exp(-I*4.*1.27*L2*length[z])+conj(Ut3f)*Um3f*exp(-I*4.*1.27*L3*length[z]);
		Ptt[z] = conj(Ut1f)*Ut1f*exp(-I*4.*1.27*L1*length[z])+conj(Ut2f)*Ut2f*exp(-I*4.*1.27*L2*length[z])+conj(Ut3f)*Ut3f*exp(-I*4.*1.27*L3*length[z]);

	}


	P[0][0]=SQR_ABS(Pee[0]);
	P[0][1]=SQR_ABS(Pem[0]);
	P[0][2]=SQR_ABS(Pet[0]);

	P[1][0]=SQR_ABS(Pme[0]);
	P[1][1]=SQR_ABS(Pmm[0]);
	P[1][2]=SQR_ABS(Pmt[0]);

	P[2][0]=SQR_ABS(Pte[0]);
	P[2][1]=SQR_ABS(Ptm[0]);
	P[2][2]=SQR_ABS(Ptt[0]);

	return P[nu_i][nu_f];
	

}


double Universe::square(double x)
{
	return x*x;
}


double Universe::min(double x, double y)
{
  if (x < y)
    return x;
  else
    return y;
}

// ----------------------------------------------------------------------------
inline void zhetrd3(complex <double> A[3][3], complex <double> Q[3][3],
                    double d[3], double e[2])
// ----------------------------------------------------------------------------
// Reduces a hermitian 3x3 matrix to real tridiagonal form by applying
// (unitary) Householder transformations:
//            [ d[0]  e[0]       ]
//    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
//            [       e[1]  d[2] ]
// The function accesses only the diagonal and upper triangular parts of
// A. The access is read-only.
// ---------------------------------------------------------------------------
{
  const int n = 3;
  complex <double> u[n], q[n];
  complex <double> omega, f;
  double K, h, g;
  int i,j;
  
  // Initialize Q to the identitity matrix
#ifndef EVALS_ONLY
  for (i=0; i < n; i++)
  {
    Q[i][i] = 1.0;
    for (j=0; j < i; j++)
      Q[i][j] = Q[j][i] = 0.0;
  }
#endif

  // Bring first row and column to the desired form 
  h = SQR_ABS(A[0][1]) + SQR_ABS(A[0][2]);
  if (real(A[0][1]) > 0)
    g = -sqrt(h);
  else
    g = sqrt(h);
  e[0] = g;
  f    = g * A[0][1];
  u[1] = conj(A[0][1]) - g;
  u[2] = conj(A[0][2]);
  
  omega = h - f;
  if (real(omega) > 0.0)
  {
    omega = 0.5 * (1.0 + conj(omega)/omega) / real(omega);
    K = 0.0;
    for (i=1; i < n; i++)
    {
      f    = conj(A[1][i]) * u[1] + A[i][2] * u[2];
      q[i] = omega * f;                  // p
      K   += real(conj(u[i]) * f);      // u* A u
    }
    K *= 0.5 * SQR_ABS(omega);

    for (i=1; i < n; i++)
      q[i] = q[i] - K * u[i];
    
    d[0] = real(A[0][0]);
    d[1] = real(A[1][1]) - 2.0*real(q[1]*conj(u[1]));
    d[2] = real(A[2][2]) - 2.0*real(q[2]*conj(u[2]));
    
    // Store inverse Householder transformation in Q
#ifndef EVALS_ONLY
    for (j=1; j < n; j++)
    {
      f = omega * conj(u[j]);
      for (i=1; i < n; i++)
        Q[i][j] = Q[i][j] - f*u[i];
    }
#endif

    // Calculate updated A[1][2] and store it in f
    f = A[1][2] - q[1]*conj(u[2]) - u[1]*conj(q[2]);
  }
  else
  {
    for (i=0; i < n; i++)
      d[i] = real(A[i][i]);
    f = A[1][2];
  }

  // Make (23) element real
  e[1] = abs(f);
#ifndef EVALS_ONLY
  if (e[1] != 0.0)
  {
    f = conj(f) / e[1];
    for (i=1; i < n; i++)
      Q[i][n-1] = Q[i][n-1] * f;
  }
#endif
}



/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


// ----------------------------------------------------------------------------
int zheevq3(complex <double> A[3][3], complex <double> Q[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a hermitian 3x3
// matrix A using the QL algorithm with implicit shifts, preceded by a
// Householder reduction to real tridiagonal form.
// The function accesses only the diagonal and upper triangular parts of A.
// The access is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The hermitian input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
// Dependencies:
//   zhetrd3()
// ----------------------------------------------------------------------------
{
  const int n = 3;
  double e[3];                 // The third element is used only as temporary workspace
  double g, r, p, f, b, s, c;  // Intermediate storage
  complex <double> t;
  int nIter;
  int m;
  int l,i,j,k;

  // Transform A to real tridiagonal form by the Householder method
  zhetrd3(A, Q, w, e);
  
  // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
  // with the QL method
  //
  // Loop over all off-diagonal elements
  for (l=0; l < n-1; l++)
  {
    nIter = 0;
    while (1)
    {
      // Check for convergence and exit iteration loop if off-diagonal
      // element e(l) is zero
      for (m=l; m <= n-2; m++)
      {
        g = fabs(w[m])+fabs(w[m+1]);
        if (fabs(e[m]) + g == g)
          break;
      }
      if (m == l)
        break;
      
      if (nIter++ >= 30)
        return -1;

      // Calculate g = d_m - k
      g = (w[l+1] - w[l]) / (e[l] + e[l]);
      r = sqrt(SQR(g) + 1.0);
      if (g > 0)
        g = w[m] - w[l] + e[l]/(g + r);
      else
        g = w[m] - w[l] + e[l]/(g - r);

      s = c = 1.0;
      p = 0.0;
      for (i=m-1; i >= l; i--)
      {
        f = s * e[i];
        b = c * e[i];
        if (fabs(f) > fabs(g))
        {
          c      = g / f;
          r      = sqrt(SQR(c) + 1.0);
          e[i+1] = f * r;
          c     *= (s = 1.0/r);
        }
        else
        {
          s      = f / g;
          r      = sqrt(SQR(s) + 1.0);
          e[i+1] = g * r;
          s     *= (c = 1.0/r);
        }
        
        g = w[i+1] - p;
        r = (w[i] - g)*s + 2.0*c*b;
        p = s * r;
        w[i+1] = g + p;
        g = c*r - b;

        // Form eigenvectors
#ifndef EVALS_ONLY
        for (k=0; k < n; k++)
        {
          t = Q[k][i+1];
          Q[k][i+1] = s*Q[k][i] + c*t;
          Q[k][i]   = c*Q[k][i] - s*t;
        }
#endif 
      }
      w[l] -= p;
      e[l]  = g;
      e[m]  = 0.0;
    }
  }

  return 0;
}





















