/* Reference: John W. RUDNICKI : FLUID MASS SOURCES AND POINT FORCES IN LINEAR ELASTIC DIFFUSIVE SOLIDS */

#include "InitialStressXYPressureNetworkLocal.h"
#include <cmath>

#define EULER 0.5772156649	/* Euler's constant gamma */
#define MAXIT 1000		/* Maximum allowed Real of iterations. */
#define FPMIN 1.0e-30	/* close to the smallest representable floting-point Real. */
#define EPS 6.0e-8		/* Desired relative error, not smaller than the machine precision. */

registerMooseObject("farmsApp", InitialStressXYPressureNetworkLocal);

InputParameters
InitialStressXYPressureNetworkLocal::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>(          "flux_q", "flux");
  params.addRequiredParam<Real>(   "density_rho_0", "fluid density");
  params.addRequiredParam<Real>(  "permeability_k", "permeability");
  params.addRequiredParam<Real>(   "viscosity_eta", "viscosity");
  params.addRequiredParam<Real>( "biotcoeff_alpha", "biot coefficient");
  params.addRequiredParam<Real>(  "undrained_nu_u", "undrained poisson's ratio");
  params.addRequiredParam<Real>("shear_modulus_mu", "shear modulus");
  params.addRequiredParam<Real>(      "drained_nu", "drained poisson's ratio");
  params.addParam<Real>("tini","initial time");
  return params;
}

InitialStressXYPressureNetworkLocal::InitialStressXYPressureNetworkLocal(const InputParameters & parameters)
  : Function(parameters),
  _flux_q(getParam<Real>("flux_q")),
  _density_rho_0(getParam<Real>("density_rho_0")),
  _permeability_k(getParam<Real>("permeability_k")),
  _viscosity_eta(getParam<Real>("viscosity_eta")),
  _biotcoeff_alpha(getParam<Real>("biotcoeff_alpha")),
  _undrained_nu_u(getParam<Real>("undrained_nu_u")),
  _shear_modulus_mu(getParam<Real>("shear_modulus_mu")),
  _drained_nu(getParam<Real>("drained_nu")),
  _tini(getParam<Real>("tini"))
{
}

/*********************************************************************
   Returns the exponential integral function
   E_n(x) = int_1^infinity e^{-x*t}/t^n dt,     for x > 0.
   C.A. Bertulani        May/15/2000
*********************************************************************/
Real expint_networklocal(int n, Real x)
{
	int i,ii,nm1;
	Real a,b,c,d,del,fact,h,psi,ans;

	nm1=n-1;
	if (n < 0 || x < 0.0 || (x==0.0 && (n==0 || n==1)))
	//err << "\n Bad arguments in expint";
	return 120e6 * 0.9;
  else {
		if (n == 0) ans=exp(-x)/x;   /* Special case */
		else {
			if (x == 0.0) ans=1.0/nm1;  /* Another special case */

			else {
				if (x > 1.0) {		/* Lentz's algorithm */
					b=x+n;
					c=1.0/FPMIN;
					d=1.0/b;
					h=d;
					for (i=1;i<=MAXIT;i++) {
						a = -i*(nm1+i);
						b += 2.0;
						d=1.0/(a*d+b);	/* Denominators cannot be zero */
						c=b+a/c;
						del=c*d;
						h *= del;
						if (fabs(del-1.0) < EPS) {
							ans=h*exp(-x);
							return ans;
						}
					}
					err << "\n Continued fraction failed in expint";
				} else {
					ans = (nm1!=0 ? 1.0/nm1 : -log(x)-EULER);	/* Set first term */
					fact=1.0;
					for (i=1;i<=MAXIT;i++) {
						fact *= -x/i;
						if (i != nm1) del = -fact/(i-nm1);
						else {
							psi = -EULER;  /* Compute psi(n) */
							for (ii=1;ii<=nm1;ii++) psi += 1.0/ii;
							del=fact*(-log(x)+psi);
						}
						ans += del;
						if (fabs(del) < fabs(ans)*EPS) return ans;
					}
					err << "\n series failed in expint";
				}
			}
		}
	}
	return ans;
}

/************************************************************************
   Returns the exponential integral function
   E_i(x) = - int_x^infinity e^{-t}/t dt = int_(-infinity)^x e^{-t}/t dt,     
   for x > 0.
   C.A. Bertulani        May/15/2000
************************************************************************/
Real ei_networklocal(Real x)
{
	int k;
	Real fact,prev,sum,term;

	if (x <= 0.0) err << "\n Bad argument in ei";
	if (x < FPMIN) return log(x)+EULER;	/* Special case: avoid failure of convergence */
	if (x <= -log(EPS)) {				/* test because of underflow.  */ 
		sum=0.0;				/* Use poer series  */
		fact=1.0;
		for (k=1;k<=MAXIT;k++) {
			fact *= x/k;
			term=fact/k;
			sum += term;
			if (term < EPS*sum) break;
		}
		if (k > MAXIT) err << "\n Series failed in ei";
		return sum+log(x)+EULER;
	} else {			/* Use asymptotic series. */
		sum=0.0;		/* Start with second term. */
		term=1.0;
		for (k=1;k<=MAXIT;k++) {
			prev=term;
			term *= k/x;
			if (term < EPS) break;
			/* Since final sum is greater than one, term itself approximates the */
			/* relative error. */
			if (term < prev) sum += term;		/* Still converging: add new term. */
			else {
				sum -= prev;		/* Diverging: subtract previous term and exit. */
				break;
			}
		}
		return exp(x)*(1.0+sum)/x;
	}
}

Real
InitialStressXYPressureNetworkLocal::value(Real t, const Point & p) const
{

  //Parameters
  //Define pi
  Real pi = 3.14159265358979323846;

  //compute R
  Real x_center = 235;
  Real y_center = 94;
  Real x_coord = p(0) - x_center; //along the strike direction
  Real y_coord = p(1) - y_center; //along the normal direction
  Real R = sqrt(x_coord*x_coord+y_coord*y_coord); //assume injection location is (0,0)

  //initialize pressure
  Real pressure = 0.0;

  //undrained lame constant
  Real drained_lambda   = 2 * _shear_modulus_mu *     _drained_nu / ( 1 - 2 *     _drained_nu );
  Real undrained_lambda = 2 * _shear_modulus_mu * _undrained_nu_u / ( 1 - 2 * _undrained_nu_u );

  //hydraulic diffusivity
  Real c = ( _permeability_k * ( undrained_lambda - drained_lambda ) * ( drained_lambda + 2 * _shear_modulus_mu ) ) / ( _viscosity_eta * _biotcoeff_alpha * _biotcoeff_alpha * ( undrained_lambda + 2 * _shear_modulus_mu ) );

  //Define z
  //Real z = R * R / ( 4 * c * (t + _tini) );

  //Compute exp integral
  //Real expIntz = expint_networklocal(1, z);

  //compute pressure
  //pressure = ( _flux_q * _viscosity_eta ) / ( 4 * pi * _density_rho_0 * _permeability_k ) * expIntz;

  //compute pressure
  pressure = ( _flux_q * _viscosity_eta ) / ( 4 * pi * _density_rho_0 * _permeability_k * R ) * std::erfc(0.5 * R / std::sqrt(c*(t+_tini)));
  
  //cap pressure
  if ( pressure > 30e6 ){
	pressure = 30e6;
  }

  return pressure;

}

#undef EPS
#undef EULER
#undef MAXIT
#undef FPMIN