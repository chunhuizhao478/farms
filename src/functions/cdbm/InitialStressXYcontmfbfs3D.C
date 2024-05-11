#include "InitialStressXYcontmfbfs3D.h"

registerMooseObject("farmsApp", InitialStressXYcontmfbfs3D);

InputParameters
InitialStressXYcontmfbfs3D::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>("maximum_value", "Maximum stress value at the bottom of simulation domain along dip");
  params.addRequiredParam<Real>("length_z", "domain length along dip");
  params.addRequiredParam<Real>("nucl_loc_x", "nucleation location x coordinate");
  params.addRequiredParam<Real>("nucl_loc_y", "nucleation location y coordinate");
  params.addRequiredParam<Real>("nucl_loc_z", "nucleation location z coordinate");
  params.addRequiredParam<Real>("nucl_patch_size", "nucleation location patch size");
  return params;
}

InitialStressXYcontmfbfs3D::InitialStressXYcontmfbfs3D(const InputParameters & parameters)
  : Function(parameters),
  _maximum_value(getParam<Real>("maximum_value")),
  _length_z(getParam<Real>("length_z")),
  _nucl_loc_x(getParam<Real>("nucl_loc_x")),
  _nucl_loc_y(getParam<Real>("nucl_loc_y")),
  _nucl_loc_z(getParam<Real>("nucl_loc_z")),
  _nucl_patch_size(getParam<Real>("nucl_patch_size")) 
{
}

Real
InitialStressXYcontmfbfs3D::value(Real /*t*/, const Point & p) const
{
  
  Real x_coord = p(0);
  Real y_coord = p(1);
  Real z_coord = p(2); //along the slip direction
  Real To = 0;

  if ( (x_coord >= _nucl_loc_x - _nucl_patch_size / 2 ) && (x_coord <= _nucl_loc_x + _nucl_patch_size / 2 ) && (y_coord >= _nucl_loc_y - _nucl_patch_size / 2 ) && (y_coord <= _nucl_loc_y + _nucl_patch_size / 2 ) && (z_coord >= _nucl_loc_z - _nucl_patch_size / 2 ) && (z_coord <= _nucl_loc_z + _nucl_patch_size / 2 ) )
  {
    //shear strength = 0.677 * 60MPa = 40.62MPa
    //overstress by 1%: 40.62 * 1.05 
    To = 42.651e6;
  }
  else{

    To = _maximum_value / _length_z * std::abs(z_coord);
  
  }
  
  return To;

}