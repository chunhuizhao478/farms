/*
Define Function for Initial Shear Stress along Strike Direction
Problem-Specific: TPV101-2D
*/

#include "InitialStrikeShearStressPerturbRSF2D.h"

registerMooseObject("farmsApp", InitialStrikeShearStressPerturbRSF2D);

InputParameters
InitialStrikeShearStressPerturbRSF2D::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStrikeShearStressPerturbRSF2D::InitialStrikeShearStressPerturbRSF2D(const InputParameters & parameters)
  : Function(parameters)
{
}

double Func_F(Real r, Real R)
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

double Func_G(Real t, Real T)
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
InitialStrikeShearStressPerturbRSF2D::value(Real t, const Point & p) const
{

  Real x_coord = p(0); //along the strike direction
  Real z_coord = p(2); //along the normal direction

  //Parameter
  Real T1_perturb = 25e6; //pertubation (Pa)
  Real R = 3000; //(m)
  Real T = 1.0; //(s)
  
  //Compute radius r
  //Real x_o = 0.0;
  //Real y_o = 0.0; //2D //along dip direction (unused)

  Real r = sqrt(x_coord*x_coord+z_coord*z_coord); //2D

  //Evalute Spatial and Temporal Function
  double Val_F = Func_F(r, R);
  double Val_G = Func_G(t, T);

  //Obtain initial shear stress
  Real T1_o = T1_perturb * Val_F * Val_G;
  
  return T1_o;

}