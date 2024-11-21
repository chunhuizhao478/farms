//Implement Lysmer damper ( same as NonreflectingBC | MASTODON ) but in explicit scheme ( follow DashpotBC | MOOSE )
//Created by Amr Ibrahim, April 18, 2024

#include "GeneralNonReflectDashpotBC_y.h"

registerMooseObject("farmsApp", GeneralNonReflectDashpotBC_y);

InputParameters
GeneralNonReflectDashpotBC_y::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredParam<unsigned int>(
      "component", "The displacement component corresponding the variable this BC acts on.");
  params.addCoupledVar("disp_x", "Displacement in the y direction");
  params.addCoupledVar("disp_z", "Displacement in the Z direction");
  params.addRequiredRangeCheckedParam<Real>(
      "p_wave_speed", "p_wave_speed>0.0", "P-wave speed of the material.");
  params.addRequiredRangeCheckedParam<Real>(
      "shear_wave_speed", "shear_wave_speed>0.0", "shear wave speed of the material.");

  return params;
}

GeneralNonReflectDashpotBC_y::GeneralNonReflectDashpotBC_y(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _disp_x_dot(isCoupled("disp_x") ? coupledDot("disp_x") : _zero),
    _d_disp_x_dot(isCoupled("disp_x") ? coupledDotDu("disp_x") : _zero),
    _disp_x_id(coupled("disp_x")),
    _disp_z_dot(isCoupled("disp_z") ? coupledDot("disp_z") : _zero),
    _d_disp_z_dot(isCoupled("disp_z") ? coupledDotDu("disp_z") : _zero),
    _disp_z_id(coupled("disp_z")),
    _u_dot(_var.uDot()), 
    _du_dot_du(_var.duDotDu()),
    _density(getMaterialPropertyByName<Real>("density")),
    _p_wave_speed(getParam<Real>("p_wave_speed")),
    _shear_wave_speed(getParam<Real>("shear_wave_speed"))
{
}

Real
GeneralNonReflectDashpotBC_y::computeQpResidual()
{
  std::vector<Real> velocity(3, 0.0);
  velocity[0] = _disp_x_dot[_qp];
  velocity[1] = _u_dot[_qp];
  velocity[2] = _disp_z_dot[_qp];

  Real normal_vel = 0.0;
  for (unsigned int i = 0; i < 3; ++i){
    normal_vel += velocity[i] * _normals[_qp](i);
  }

  return _test[_i][_qp] * _density[_qp] *
         (_p_wave_speed * normal_vel * _normals[_qp](_component) +
          _shear_wave_speed * (velocity[_component] - normal_vel * _normals[_qp](_component)));
}
Real
GeneralNonReflectDashpotBC_y::computeQpJacobian()
{

  return _test[_i][_qp] * _density[_qp] * _phi[_j][_qp] *
         (_p_wave_speed * _du_dot_du[_qp] * _normals[_qp](_component) +
          _shear_wave_speed * (_du_dot_du[_qp] - _du_dot_du[_qp] * _normals[_qp](_component)));
}
Real
GeneralNonReflectDashpotBC_y::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _disp_x_id)
    return _test[_i][_qp] * _density[_qp] *
         (_p_wave_speed *  _d_disp_x_dot[_qp] * _normals[_qp](_disp_x_id) * _normals[_qp](_component) +
          _shear_wave_speed * ( _d_disp_x_dot[_qp] -  _d_disp_x_dot[_qp] * _normals[_qp](_disp_x_id) * _normals[_qp](_component)))* _phi[_j][_qp];
  else if (jvar == _disp_z_id)
   return _test[_i][_qp] * _density[_qp] *
         (_p_wave_speed *  _d_disp_z_dot[_qp] * _normals[_qp](_disp_z_id) * _normals[_qp](_component) +
          _shear_wave_speed * ( _d_disp_z_dot[_qp] -  _d_disp_z_dot[_qp] * _normals[_qp](_disp_z_id) * _normals[_qp](_component)))* _phi[_j][_qp];
  return 0.0;
}

