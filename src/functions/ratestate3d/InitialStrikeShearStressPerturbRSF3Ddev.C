/*
Define Function for Initial Shear Stress along Strike Direction
Problem-Specific: TPV101-3D
*/

#include "InitialStrikeShearStressPerturbRSF3Ddev.h"

registerMooseObject("farmsApp", InitialStrikeShearStressPerturbRSF3Ddev);

InputParameters
InitialStrikeShearStressPerturbRSF3Ddev::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>("stress_perturb","stress perturbation");
  return params;
}

InitialStrikeShearStressPerturbRSF3Ddev::InitialStrikeShearStressPerturbRSF3Ddev(const InputParameters & parameters)
  : Function(parameters),
  _stress_perturb(getParam<Real>("stress_perturb"))
{
}

double Function_F_dev(Real r, Real R)
{
    double Val_F = 0.0;
    if (r < R){
        Val_F = exp((r*r)/(r*r-R*R));
    }
    else{
        Val_F = 0.0;
    }
    return Val_F;
}

double Function_G_dev(Real t, Real T)
{
    double Val_G = 0.0;
    if (t < T && t > 0 ){
        Val_G = exp(((t-T)*(t-T))/(t*(t-2*T)));
    }
    else if (t >=T){
        Val_G = 1;
    }
    return Val_G;
}

Real
InitialStrikeShearStressPerturbRSF3Ddev::value(Real t, const Point & p) const
{

  Real x_coord = p(0); //along the strike direction
  Real z_coord = p(2); //along the dip direction

  //Parameter
  Real T1_perturb = _stress_perturb; //pertubation (Pa)
  Real R = 3000; //(m)
  Real T = 1.0; //(s)
  
  //Compute radius r
  Real x_o = 0.0; //3D //along strike direction
  Real z_o = 0.0; //3D //along dip direction

  Real r = sqrt((x_coord-x_o)*(x_coord-x_o)+(z_coord-z_o)*(z_coord-z_o)); //3D

  //Evalute Spatial and Temporal Function
  double Val_F = Function_F_dev(r, R);
  double Val_G = Function_G_dev(t, T);

  //Obtain initial shear stress
  Real T1_o = T1_perturb * Val_F * Val_G;
  
  return T1_o;

}