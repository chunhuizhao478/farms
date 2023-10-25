#include "InitialStressXYcontmfbfs.h"

registerMooseObject("farmsApp", InitialStressXYcontmfbfs);

InputParameters
InitialStressXYcontmfbfs::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressXYcontmfbfs::InitialStressXYcontmfbfs(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressXYcontmfbfs::value(Real /*t*/, const Point & p) const
{

  // Real T1_o = 0.0;

  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction

  // //Parameter
  // Real T1_perturb = 11.6e6; //pertubation (Pa) 81.6e6-70e6
  // Real R = 500; //(m) !!!

  // Real r = sqrt(x_coord*x_coord+y_coord*y_coord); //2D

  // Real Val_F = 0.0;
  // if (r < R){
  //   Val_F = exp((r*r)/(r*r-R*R));
  // }
  // else{
  //   Val_F = 0.0;
  // }

  // Real T1_perturb_o = T1_perturb * Val_F;

  // T1_o = 70e6 + T1_perturb_o;

  Real T1_o = 0;

  if (x_coord >= -100 and x_coord <= 100 and y_coord >= -100 and y_coord <= 100)
  {
    T1_o = 81.6e6;
  }
  else{
    T1_o = 70e6;
  }
  
  return T1_o;

}