#include "InitialStresscontmfbfs3D.h"

registerMooseObject("farmsApp", InitialStresscontmfbfs3D);

InputParameters
InitialStresscontmfbfs3D::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>("maximum_value", "Maximum stress value at the bottom of simulation domain along dip");
  params.addRequiredParam<Real>("length_z", "domain length along dip");
  return params;
}

InitialStresscontmfbfs3D::InitialStresscontmfbfs3D(const InputParameters & parameters)
  : Function(parameters),
  _maximum_value(getParam<Real>("maximum_value")),
  _length_z(getParam<Real>("length_z"))
{
}

Real
InitialStresscontmfbfs3D::value(Real /*t*/, const Point & p) const
{
  
  Real z_coord = p(2); //along the slip direction
  Real To = 0;

  //add 1e-3 small number to avoid non-zero xi value
  To = 1e-3 + _maximum_value / _length_z * std::abs(z_coord);
  
  return To;

}