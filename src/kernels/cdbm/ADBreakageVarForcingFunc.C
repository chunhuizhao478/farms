/*
Update Breakage Variable Using Runge-Kutta own solver

- 10/5/2023 Chunhui Zhao

Include power-law correction on Cd (function of strain rate e)
if (e < 1e-4){ Cd = 10 };
else{ Cd = std::std::std::pow(10, log10(1+m*log10(e/1e-4)) ) }

Cb = CdCb_multiplier * Cd

*/

#include "ADBreakageVarForcingFunc.h"

registerMooseObject("farmsApp", ADBreakageVarForcingFunc);

InputParameters
ADBreakageVarForcingFunc::validParams()
{
  InputParameters params = ADKernel::validParams();

  //constant parameters
  params.addRequiredParam<Real>(  "CBCBH_multiplier", "coefficient of healing for breakage evolution");
  params.addRequiredParam<Real>(         "a0", "parameters in granular states");
  params.addRequiredParam<Real>(         "a1", "parameters in granular states");
  params.addRequiredParam<Real>(         "a2", "parameters in granular states");
  params.addRequiredParam<Real>(         "a3", "parameters in granular states");
  params.addRequiredParam<Real>(     "xi_min", "strain invariant ratio at minimum value");
  params.addRequiredParam<Real>(     "xi_max", "strain invariant ratio at maximum value");
  params.addRequiredParam<Real>(       "xi_d", "strain invariant ratio at maximum value");
  params.addRequiredParam<Real>(       "xi_0", "strain invariants ratio: onset of damage evolution");
  params.addRequiredParam<Real>(       "xi_1", "strain invariants ratio: onset of breakage healing");
  params.addRequiredParam<Real>( "beta_width", "coefficient gives width of transitional region");
  params.addRequiredParam<Real>("CdCb_multiplier", "multiplier between Cd and Cb");

  //variable parameters
  params.addRequiredCoupledVar(  "alpha_old", "damage variable at previous time step");
  params.addRequiredCoupledVar(      "B_old", "breakage variable at previous time step");
  params.addRequiredCoupledVar(     "xi_old", "strain invariant ratio at previous time step");
  params.addRequiredCoupledVar(     "I2_old", "second strain invariant at previous time step");
  params.addRequiredCoupledVar(     "mu_old", "shear modulus at previous time step");
  params.addRequiredCoupledVar( "lambda_old", "first lame constant at previous time step");
  params.addRequiredCoupledVar(  "gamma_old", "damage modulus at previous time step");

  params.addParam<Real>( "Cd_constant", 0.0, "constant Cd value for option 2 only");

  return params;
}

ADBreakageVarForcingFunc::ADBreakageVarForcingFunc(const InputParameters & parameters)
 : ADKernel(parameters),
  _CBCBH_multiplier(getParam<Real>("CBCBH_multiplier")),
  _a0(getParam<Real>("a0")),
  _a1(getParam<Real>("a1")),
  _a2(getParam<Real>("a2")),
  _a3(getParam<Real>("a3")),
  _xi_min(getParam<Real>("xi_min")),
  _xi_max(getParam<Real>("xi_max")),
  _xi_d(getParam<Real>("xi_d")),
  _xi_0(getParam<Real>("xi_0")),
  _xi_1(getParam<Real>("xi_1")),
  _beta_width(getParam<Real>("beta_width")),
  _CdCb_multiplier(getParam<Real>("CdCb_multiplier")),
  _alpha_old(adCoupledValue("alpha_old")),
  _B_old(adCoupledValue("B_old")),
  _xi_old(adCoupledValue("xi_old")),
  _I2_old(adCoupledValue("I2_old")),
  _mu_old(adCoupledValue("mu_old")),
  _lambda_old(adCoupledValue("lambda_old")),
  _gamma_old(adCoupledValue("gamma_old")),
  _Cd_constant(getParam<Real>("Cd_constant"))
{
}

ADReal
ADBreakageVarForcingFunc::computeQpResidual()
{ 
  
    //get parameters
    ADReal alpha = _alpha_old[_qp];
    ADReal B = _B_old[_qp];
    ADReal I2 = _I2_old[_qp];
    ADReal xi = _xi_old[_qp];

    //Initialize Cd
    ADReal Cd = _Cd_constant;

    //Compute C_B
    ADReal C_B = _CdCb_multiplier * Cd;

    //Compute C_BH //keep constant 
    ADReal C_BH = _CBCBH_multiplier * C_B;

    //
    ADReal alphacr = computeAlphaCr(xi);
    ADReal Prob = 1.0 / ( std::exp( (alphacr - alpha) / _beta_width ) + 1.0 );

    //no healing //this formulation is used in the splitstrain article
    if ( xi >= _xi_0 && xi <= _xi_max ){
        return -1.0 * C_B * Prob * (1-B) * I2 * (xi - _xi_0) * _test[_i][_qp]; //could heal if xi < xi_0
    }
    else if ( xi < _xi_0 && xi >= _xi_min ){
        return -1.0 * C_BH * I2 * ( xi - _xi_0 ) * _test[_i][_qp];
    }
    else{
        mooseError("xi_old is OUT-OF-RANGE!.");
        return 0;
    }

}

/// Function: Compute alpha_cr based on the current xi
ADReal 
ADBreakageVarForcingFunc::computeAlphaCr(ADReal xi)
{
  ADReal alphacr;
  if ( xi < _xi_0 )
  {
    alphacr = 1.0;
  } 
  else if ( xi > _xi_0 && xi <= _xi_1 )
  {
    alphacr = ((xi*2.76e5-7.100521107637101e2*xi*7.5e2-7.100521107637101e2*1.4e3-7.100521107637101e+2*std::pow(xi,3)*1.25e2+std::pow(xi,3)*4.6e4+std::sqrt((7.100521107637101e2*3.68e2-3.19799e5)*(xi*(-1.44e3)-std::pow(xi,2)*2.1e3+std::pow(xi,3)*5.6e2+std::pow(xi,4)*3.0e2+std::pow(xi,6)*2.5e1+3.576e3)*(-3.590922148807814e-1))*5.9e1+5.152e5)*(-5.9e1/4.0))/(xi*3.837588e7-7.100521107637101e2*xi*4.416e4+7.100521107637101e2*4.048e3-7.100521107637101e2*std::pow(xi,2)*2.76e4+std::pow(xi,2)*2.3984925e7-3.517789e6);
  }
  else if ( xi > _xi_1 && xi <= _xi_max )
  {
    alphacr = 6.408e10/(7.100521107637101e2*1.737762711864407e8+xi*(7.100521107637101e2*1.086101694915254e8-3.996854237288136e10)-6.394966779661017e10);
  }
  else
  {
    std::cout<<"xi: "<<xi<<std::endl;
    mooseError("xi exceeds the maximum allowable range!");
  }
  return alphacr;

}