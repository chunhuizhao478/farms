/*
Define Function for Initial Shear Stress along Strike Direction
Problem-Specific: TPV101-2D
*/

#include "InitialStrikeShearStressPerturbRSF2Ddev.h"

registerMooseObject("farmsApp", InitialStrikeShearStressPerturbRSF2Ddev);

InputParameters
InitialStrikeShearStressPerturbRSF2Ddev::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>("stress_perturb","stress perturbation value");
  return params;
}

InitialStrikeShearStressPerturbRSF2Ddev::InitialStrikeShearStressPerturbRSF2Ddev(const InputParameters & parameters)
  : Function(parameters),
  _stress_perturb(getParam<Real>("stress_perturb"))
{
}

Real Func_F_dev(Real r, Real R)
{
    Real Val_F = 0.0;
    if (r < R){
        Val_F = exp((r*r)/(r*r-R*R));
    }
    else{
        Val_F = 0.0;
    }
    return Val_F;
}

Real Func_G_dev(Real t, Real T)
{
    Real Val_G = 0.0;
    if (t < T && t > 0 ){
        Val_G = exp(((t-T)*(t-T))/(t*(t-2*T)));
    }
    else if (t >=T){
        Val_G = 1;
    }
    return Val_G;
}

Real
InitialStrikeShearStressPerturbRSF2Ddev::value(Real t, const Point & p) const
{

  Real x_coord = p(0); //along the strike direction
  Real z_coord = p(2); //along the normal direction

  //Parameter
  Real T1_perturb = _stress_perturb; //pertubation (Pa)
  Real R = 3000; //(m)
  Real T = 1.0; //(s)
  
  //Compute radius r
  //Real x_o = 0.0;
  //Real y_o = 0.0; //2D //along dip direction (unused)

  Real r = sqrt(x_coord*x_coord+z_coord*z_coord); //2D

  //Evalute Spatial and Temporal Function
  Real Val_F = Func_F_dev(r, R);
  Real Val_G = Func_G_dev(t, T);

  //Obtain initial shear stress
  Real T1_o = T1_perturb * Val_F * Val_G;
  
  return T1_o;

}