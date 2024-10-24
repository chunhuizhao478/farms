#include "BreakageVarUpdate.h"

registerMooseObject("farmsApp", BreakageVarUpdate);

InputParameters
BreakageVarUpdate::validParams()
{
  InputParameters params = AuxKernel::validParams();

  //constant parameters
  params.addParam<Real>(        "C_B", 1.0, "coefficient gives positive breakage evolution");
  params.addParam<Real>(       "C_BH", 1.0, "coefficient of healing for breakage evolution");
  params.addParam<Real>(         "a0", 1.0, "parameters in granular states");
  params.addParam<Real>(         "a1", 1.0, "parameters in granular states");
  params.addParam<Real>(         "a2", 1.0, "parameters in granular states");
  params.addParam<Real>(         "a3", 1.0, "parameters in granular states");
  params.addParam<Real>(     "xi_min", 1.0, "strain invariant ratio at minimum value");
  params.addParam<Real>(     "xi_max", 1.0, "strain invariant ratio at maximum value");
  params.addParam<Real>(       "xi_d", 1.0, "strain invariant ratio at maximum value");
  params.addParam<Real>(       "xi_0", 1.0, "strain invariants ratio: onset of damage evolution");
  params.addParam<Real>(       "xi_1", 1.0, "strain invariants ratio: onset of breakage healing");
  params.addParam<Real>( "beta_width", 1.0, "coefficient gives width of transitional region");

  //variable parameters
  params.addRequiredCoupledVar(  "alpha_old", "damage variable at previous time step");
  params.addRequiredCoupledVar(      "B_old", "breakage variable at previous time step");
  params.addRequiredCoupledVar(     "xi_old", "strain invariant ratio at previous time step");
  params.addRequiredCoupledVar(     "I2_old", "second strain invariant at previous time step");
  params.addRequiredCoupledVar(     "mu_old", "shear modulus at previous time step");
  params.addRequiredCoupledVar( "lambda_old", "first lame constant at previous time step");
  params.addRequiredCoupledVar(  "gamma_old", "damage modulus at previous time step");

  return params;
}

BreakageVarUpdate::BreakageVarUpdate(const InputParameters & parameters)
  : AuxKernel(parameters),
  _C_B(getParam<Real>("C_B")),
  _C_BH(getParam<Real>("C_BH")),
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
  _alpha_old(coupledValue("alpha_old")),
  _B_old(coupledValue("B_old")),
  _xi_old(coupledValue("xi_old")),
  _I2_old(coupledValue("I2_old")),
  _mu_old(coupledValue("mu_old")),
  _lambda_old(coupledValue("lambda_old")),
  _gamma_old(coupledValue("gamma_old"))
{
}

Real
BreakageVarUpdate::computeValue()
{   
    //get parameters
    Real alpha = _alpha_old[_qp];
    Real B = _B_old[_qp];
    Real I2 = _I2_old[_qp];
    Real xi = _xi_old[_qp];
    Real mu = _mu_old[_qp];
    Real gamma_damaged = _gamma_old[_qp];
    Real lambda = _lambda_old[_qp];

    //ode solve
    Real B_update = OdeIntegrator(alpha,B,I2,xi,mu,gamma_damaged,lambda);

    return B_update;
}

///Helper Function
///Function: Dormand-Prince Method ode integrator
Real 
BreakageVarUpdate::OdeIntegrator(Real alpha,
                                 Real B,
                                 Real I2,
                                 Real xi,
                                 Real mu,
                                 Real gamma_damaged,
                                 Real lambda)
{
  Real ynew = 0;
  // //initialize k
  // Real k1 = 0.0; Real k2 = 0.0; Real k3 = 0.0;
  // Real k4 = 0.0; Real k5 = 0.0; Real k6 = 0.0; Real k7 = 0.0;
  // //parameters
  // Real a1 = 1.0/5.0;
  // Real b1 = 3.0/40.0;       Real b2 = 9.0/40.0; 
  // Real c1 = 44.0/45.0;      Real c2 = -56.0/15.0;      Real c3 = 32.0/9.0; 
  // Real d1 = 19372.0/6561.0; Real d2 = -25360.0/2187.0; Real d3 = 64448.0/6561.0;  Real d4 = -212.0/729.0; 
  // Real e1 = 9017.0/3168.0;  Real e2 = -355.0/33.0;     Real e3 = 46732.0/5247.0;  Real e4 = 49.0/176.0;     Real e5 = -5103.0/18656.0; 
  // Real f1 = 35.0/384.0;     Real f3 = 500.0/1113.0;    Real f4 = 125.0/192.0;     Real f5 = -2187.0/6784.0; Real f6 = 11.0/84.0;
  // //
  // Real y1 = 35.0/384.0; Real y3 = 500.0/1113.0; Real y4 = 125.0/192.0;
  // Real y5 = -2187.0/6784.0; Real y6 = 11.0/84.0;

  // k1 = _dt * computeBreakageEvolution(alpha                                                  ,B                                                  ,I2,xi,mu,gamma_damaged,lambda);
  // k2 = _dt * computeBreakageEvolution(alpha + a1 * k1                                        ,B + a1 * k1                                        ,I2,xi,mu,gamma_damaged,lambda);
  // k3 = _dt * computeBreakageEvolution(alpha + b1 * k1 + b2 * k2                              ,B + b1 * k1 + b2 * k2                              ,I2,xi,mu,gamma_damaged,lambda);
  // k4 = _dt * computeBreakageEvolution(alpha + c1 * k1 + c2 * k2 + c3 * k3                    ,B + c1 * k1 + c2 * k2 + c3 * k3                    ,I2,xi,mu,gamma_damaged,lambda);
  // k5 = _dt * computeBreakageEvolution(alpha + d1 * k1 + d2 * k2 + d3 * k3 + d4 * k4          ,B + d1 * k1 + d2 * k2 + d3 * k3 + d4 * k4          ,I2,xi,mu,gamma_damaged,lambda);
  // k6 = _dt * computeBreakageEvolution(alpha + e1 * k1 + e2 * k2 + e3 * k3 + e4 * k4 + e5 * k5,B + e1 * k1 + e2 * k2 + e3 * k3 + e4 * k4 + e5 * k5,I2,xi,mu,gamma_damaged,lambda);
  // k7 = _dt * computeBreakageEvolution(alpha + f1 * k1 + f3 * k3 + f4 * k4 + f5 * k5 + f6 * k6,B + f1 * k1 + f3 * k3 + f4 * k4 + f5 * k5 + f6 * k6,I2,xi,mu,gamma_damaged,lambda);
  
  //ynew = B + y1 * k1 + y3 * k3 + y4 * k4 + y5 * k5 + y6 * k6; 

  Real k1 = 0.0;
  k1 = _dt * computeBreakageEvolution(alpha,B,I2,xi,mu,gamma_damaged,lambda);
  ynew = B + k1; //use first-order approx

  return ynew;
}

Real
BreakageVarUpdate::computeBreakageEvolution(Real alpha,
                                            Real B,
                                            Real I2,
                                            Real xi,
                                            Real /*mu*/,
                                            Real /*gamma_damaged*/,
                                            Real /*lambda*/)
{ 
  Real Prob = computeGranularStateProb(alpha, xi);

  // //with healing
  // if ( xi >= _xi_d && xi <= _xi_max ){
  //   return _C_B * Prob * (1-B) * I2 * ( ( mu - gamma_damaged * xi + 0.5 * lambda * pow(xi,2) ) -  ( _a0 + _a1 * xi + _a2 * pow(xi,2) + _a3 * pow(xi,3) ) );
  // }
  // else if ( xi < _xi_d && xi >= _xi_min ){
  //   //with healing
  //   //return _C_BH * I2 * ( ( mu - gamma_damaged * xi + 0.5 * lambda * pow(xi,2) ) -  ( _a0 + _a1 * xi + _a2 * pow(xi,2) + _a3 * pow(xi,3) ) );
  //   //no healing
  //   return 0.0;
  // }
  // else{
  //   std::cout<<"xi: "<<xi<<std::endl;
  //   mooseError("xi_old is OUT-OF-RANGE!.");
  //   return 0;
  // }

  //no healing
  if ( xi >= _xi_0 ){
    return _C_B * Prob * (1-B) * I2 * (xi - _xi_0); //could heal if xi < xi_0
    //return _C_B * Prob * (1-B) * I2 * abs(xi - _xi_0);  //forbid any healing
  }
  else{
    return 0;
  }
}

/// Function: Probability for material to be in granular state
Real 
BreakageVarUpdate::computeGranularStateProb(Real alpha,
                                            Real xi)
{
  Real Prob = 0;
  Real alphacr = computeAlphaCr(xi);
  Prob = 1.0 / ( exp( (alphacr - alpha) / _beta_width ) + 1.0 );
  return Prob;
}

/// Function: Compute alpha_cr based on the current xi
Real 
BreakageVarUpdate::computeAlphaCr(Real xi)
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
    alphacr = ((xi*2.76e5-7.100521107637101e2*xi*7.5e2-7.100521107637101e2*1.4e3-7.100521107637101e+2*pow(xi,3)*1.25e2+pow(xi,3)*4.6e4+sqrt((7.100521107637101e2*3.68e2-3.19799e5)*(xi*(-1.44e3)-pow(xi,2)*2.1e3+pow(xi,3)*5.6e2+pow(xi,4)*3.0e2+pow(xi,6)*2.5e1+3.576e3)*(-3.590922148807814e-1))*5.9e1+5.152e5)*(-5.9e1/4.0))/(xi*3.837588e7-7.100521107637101e2*xi*4.416e4+7.100521107637101e2*4.048e3-7.100521107637101e2*pow(xi,2)*2.76e4+pow(xi,2)*2.3984925e7-3.517789e6);
  }
  else if ( xi > _xi_1 && xi <= _xi_max )
  {
    //lambda_o,shear_modulus_o = 1e11
    //alphacr = (1.078052110763710e+03)/(125*(5*xi + 8));
    //lambda_o,shear_modulus_o = 32.04e9
    alphacr = 6.408e10/(7.100521107637101e2*1.737762711864407e8+xi*(7.100521107637101e2*1.086101694915254e8-3.996854237288136e10)-6.394966779661017e10);
  }
  else
  {
    std::cout<<"xi: "<<xi<<std::endl;
    mooseError("xi exceeds the maximum allowable range!");
  }
  return alphacr;

}

