/*
Update Breakage Variable Using Runge-Kutta own solver

- 10/5/2023 Chunhui Zhao

Include power-law correction on Cd (function of strain rate e)
if (e < 1e-4){ Cd = 10 };
else{ Cd = pow(10, log10(1+m*log10(e/1e-4)) ) }

Cb = CdCb_multiplier * Cd

*/

#include "BreakageVarForcingFuncDev.h"

registerMooseObject("farmsApp", BreakageVarForcingFuncDev);

InputParameters
BreakageVarForcingFuncDev::validParams()
{
  InputParameters params = Kernel::validParams();

  //constant parameters
  params.addParam<Real>(   "C_d_min", 10.0, "coefficient gives positive damage evolution (small strain e < 1e-4 threshold value)");
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
  params.addRequiredParam<Real>(     "m", "Cd power-law correction index");
  params.addParam<Real>("mechanical_strain_rate_threshold", 0, "threshold value for strain rate such that Cd takes constant value Cd_min if strain rate below this value.");
  params.addRequiredParam<Real>("CdCb_multiplier", "multiplier between Cd and Cb");
  params.addParam<Real>( "scale", 1.0, "scale the Cd power-law");

  //variable parameters
  params.addRequiredCoupledVar(  "alpha_old", "damage variable at previous time step");
  params.addRequiredCoupledVar(      "B_old", "breakage variable at previous time step");
  params.addRequiredCoupledVar(     "xi_old", "strain invariant ratio at previous time step");
  params.addRequiredCoupledVar(     "I2_old", "second strain invariant at previous time step");
  params.addRequiredCoupledVar(     "mu_old", "shear modulus at previous time step");
  params.addRequiredCoupledVar( "lambda_old", "first lame constant at previous time step");
  params.addRequiredCoupledVar(  "gamma_old", "damage modulus at previous time step");
  params.addCoupledVar("mechanical_strain_rate", 0.0, "strain rate");

  //add options
  params.addRequiredParam<int>( "option", "option 1 : Cd power-law; option 2 : use constant Cd");
  params.addParam<Real>( "Cd_constant", 0.0, "constant Cd value for option 2 only");

  //add healing
  params.addParam<bool>("healing", false, "if turn on healing, true = on, false = off, default is false = off");

  return params;
}

BreakageVarForcingFuncDev::BreakageVarForcingFuncDev(const InputParameters & parameters)
 : Kernel(parameters),
  _Cd_min(getParam<Real>("C_d_min")),
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
  _scale(getParam<Real>("scale")),
  _m(getParam<Real>("m")),
  _mechanical_strain_rate_threshold(getParam<Real>("mechanical_strain_rate_threshold")),
  _CdCb_multiplier(getParam<Real>("CdCb_multiplier")),
  _alpha_old(coupledValue("alpha_old")),
  _B_old(coupledValue("B_old")),
  _xi_old(coupledValue("xi_old")),
  _I2_old(coupledValue("I2_old")),
  _mu_old(coupledValue("mu_old")),
  _lambda_old(coupledValue("lambda_old")),
  _gamma_old(coupledValue("gamma_old")),
  _mechanical_strain_rate(coupledValue("mechanical_strain_rate")),
  _option(getParam<int>("option")),
  _Cd_constant(getParam<Real>("Cd_constant")),
  _healing(getParam<bool>("healing"))
{
}

Real
BreakageVarForcingFuncDev::computeQpResidual()
{ 
  
 //get parameters
    Real alpha = _alpha_old[_qp];
    Real B = _B_old[_qp];
    Real I2 = _I2_old[_qp];
    Real xi = _xi_old[_qp];
    // Real mu = _mu_old[_qp];
    // Real gamma_damaged = _gamma_old[_qp];
    // Real lambda = _lambda_old[_qp];

    //Power-law correction
    //Initialize Cd
    Real Cd = 0;
    //Check options
    if ( _option == 1 ){

      //close _option == 1
      //mooseError("Option 1 is NOT available!");

      //power-law correction on coefficient Cd(function of strain rate)
      if ( _mechanical_strain_rate[_qp] < _mechanical_strain_rate_threshold ) //Cd follows power-law
      {
        Cd = _scale * pow(10, 1 + _m * log10( _mechanical_strain_rate[_qp] / _mechanical_strain_rate_threshold ) ) * _Cd_min;
      }
      else if ( _mechanical_strain_rate[_qp] < 0 ){ //Cd remains constant
        Cd = _Cd_min;
      }

    }
    else if ( _option == 2 ){

      if ( _Cd_constant == 0.0 ){
        mooseError("For option 2, need to provide nonzero Cd_constant value !");
      }
      else{
        Cd = _Cd_constant;
      }

    }
    else{
      mooseError("Please provide valid option number!");
    }

    //Compute C_B
    Real C_B = _CdCb_multiplier * Cd;

    //Compute C_BH //keep constant 
    Real C_BH = 1e4;

    //
    Real Prob = 0;
    Real alphacr = computeAlphaCr(xi);
    Prob = 1.0 / ( exp( (alphacr - alpha) / _beta_width ) + 1.0 );

    //no healing //this formulation is used in the splitstrain article
    if ( xi >= _xi_0 && xi <= _xi_max ){
        return -1.0 * C_B * Prob * (1-B) * I2 * (xi - _xi_0) * _test[_i][_qp]; //could heal if xi < xi_0
    }
    else if ( xi < _xi_0 && xi >= _xi_min ){
      
      if ( _healing == true ){
        return -1.0 * C_BH * I2 * ( xi - _xi_0 ) * _test[_i][_qp];
      }
      else{
        return 0.0;
      }
    }
    else{
      //mooseError("xi_old is OUT-OF-RANGE!.");
      return 0;
    }

    // else{ //with healing
    //   //add pre-check to the value (cause problem if encountered large values)
    //   if ( xi >= _xi_d && xi <= _xi_max ){
    //     return -1.0 * C_B * Prob * (1-B) * I2 * ( ( mu - gamma_damaged * xi + 0.5 * lambda * pow(xi,2) ) -  ( _a0 + _a1 * xi + _a2 * pow(xi,2) + _a3 * pow(xi,3) ) ) * _test[_i][_qp];
    //   }
    //   else if ( xi < _xi_d && xi >= _xi_min ){
    //     if ( B != 1.0 ){ //no healing until B reaches its maximum
    //       return 0.0 * _test[_i][_qp];
    //     }
    //     else{
    //       return -1.0 * C_BH * I2 * ( ( mu - gamma_damaged * xi + 0.5 * lambda * pow(xi,2) ) -  ( _a0 + _a1 * xi + _a2 * pow(xi,2) + _a3 * pow(xi,3) ) ) * _test[_i][_qp];
    //     }
    //   }
    //   else{
    //     std::cout<<"xi: "<<xi<<std::endl;
    //     mooseError("xi_old is OUT-OF-RANGE!.");
    //     return 0;
    //   }
    // }
}

Real
BreakageVarForcingFuncDev::computeQpJacobian()
{
  return 0.0;
}

/// Function: Compute alpha_cr based on the current xi
Real 
BreakageVarForcingFuncDev::computeAlphaCr(Real xi)
{
  Real alphacr;
  if ( xi < _xi_0 )
  {
    alphacr = 1.0;
  } 
  else if ( xi > _xi_0 && xi <= _xi_1 )
  {
    //lambda_o,shear_modulus_o = 1e11
    //alphacr = (2.156104221527420*(30*xi - sqrt(25*pow(xi,6) + 300*pow(xi,4) + 560*pow(xi,3) - 2100*pow(xi,2) - 1440*xi + 3576) + 5*pow(xi,3) + 56))/(75*pow(xi,2) + 120*xi - 11);
    //lambda_o,shear_modulus_o = 32.04e9
    //alphacr = ((xi*2.76e5-7.100521107637101e2*xi*7.5e2-7.100521107637101e2*1.4e3-7.100521107637101e+2*pow(xi,3)*1.25e2+pow(xi,3)*4.6e4+sqrt((7.100521107637101e2*3.68e2-3.19799e5)*(xi*(-1.44e3)-pow(xi,2)*2.1e3+pow(xi,3)*5.6e2+pow(xi,4)*3.0e2+pow(xi,6)*2.5e1+3.576e3)*(-3.590922148807814e-1))*5.9e1+5.152e5)*(-5.9e1/4.0))/(xi*3.837588e7-7.100521107637101e2*xi*4.416e4+7.100521107637101e2*4.048e3-7.100521107637101e2*pow(xi,2)*2.76e4+pow(xi,2)*2.3984925e7-3.517789e6);
    alphacr = ((xi*1.616908125e11-1.420703785875508e6*xi*1.875e5-1.420703785875508e6*4.2875e5+sqrt((1.420703785875508e6*8.62351e5-1.381024247201e12)*(xi*-4.41e4-pow(xi,2)*5.25e4+pow(xi,3)*1.715e4+pow(xi,4)*7.5e3+pow(xi,6)*6.25e2+9.6609e4)*(-1.201932276327807e-1))*5.099e3-1.420703785875508e6*pow(xi,3)*3.125e4+pow(xi,3)*2.694846875e10+3.6973299125e11)*(-2.5495e3))/(xi*5.075264108463675e15-1.420703785875508e6*xi*3.169139925e9-1.420703785875508e6*4.53596626e8-1.420703785875508e6*pow(xi,2)*1.616908125e9+pow(xi,2)*2.589420463501875e15+7.26418754027726e14);
  }
  else if ( xi > _xi_1 && xi <= _xi_max )
  {
    //lambda_o,shear_modulus_o = 1e11
    //alphacr = (1.078052110763710e+03)/(125*(5*xi + 8));
    //lambda_o,shear_modulus_o = 32.04e9
    alphacr = 6.408e10/(1.420703785875508e6*1.231582663267307e5+xi*(1.420703785875508e6*6.283585016669935e4-5.418655822710335e10)-1.062056541251226e11);
  }
  else
  {
    std::cout<<"xi: "<<xi<<std::endl;
    mooseError("xi exceeds the maximum allowable range!");
  }
  return alphacr;

}