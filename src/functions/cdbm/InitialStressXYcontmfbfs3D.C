#include "InitialStressXYcontmfbfs3D.h"

registerMooseObject("farmsApp", InitialStressXYcontmfbfs3D);

InputParameters
InitialStressXYcontmfbfs3D::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressXYcontmfbfs3D::InitialStressXYcontmfbfs3D(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressXYcontmfbfs3D::value(Real /*t*/, const Point & p) const
{

  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction
  Real z_coord = p(2); //along the slip direction

  Real T1_o = 0;

  if ((x_coord<=(0.0+100))&&(x_coord>=(0.0-100))&& (z_coord<=(-400+100))&&(z_coord>=(-400-100))&&(y_coord<=(0+100))&&(y_coord>=(0-100)))
  {
      T1_o = 81.6e6;
  }
  else{
      T1_o = 70e6;
  }
  
  return T1_o;

}