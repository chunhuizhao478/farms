#include "InitialShearStress3D.h"

registerMooseObject("farmsApp", InitialShearStress3D);

InputParameters
InitialShearStress3D::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialShearStress3D::InitialShearStress3D(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialShearStress3D::value(Real /*t*/, const Point & p) const
{
  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the dip direction
  Real z_coord = p(2); //along the normal direction

  Real T1_o = 0;

  if (x_coord >= -2870 and x_coord <= -2465 and y_coord >= -4200 and y_coord <= -3800 and z_coord >= -500 and z_coord <= 500)
  {
    T1_o = 81.6e6;
  }
  else{
    T1_o = 70e6;
  }
  
  return T1_o;

}