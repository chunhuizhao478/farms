//Implement Lysmer damper ( same as NonreflectingBC | MASTODON ) but in explicit scheme ( follow DashpotBC | MOOSE )
//Created by Amr Ibrahim, April 18, 2024

#include "BiotSlowWaveNonReflectDashpotBC_x.h"

registerMooseObject("farmsApp", BiotSlowWaveNonReflectDashpotBC_x);

InputParameters
BiotSlowWaveNonReflectDashpotBC_x::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredParam<unsigned int>(
      "component", "The Fluid velocity component corresponding the variable this BC acts on.");
  params.addCoupledVar("fluid_vel_y", "Fluid velocity in the y direction");
  params.addCoupledVar("fluid_vel_z", "Fluid velocity in the Z direction");
  params.addRequiredRangeCheckedParam<Real>(
      "_biot_p_wave_speed", "p_wave_speed>0.0", "P-wave speed of the material.");
  return params;
}

BiotSlowWaveNonReflectDashpotBC_x::BiotSlowWaveNonReflectDashpotBC_x(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _fluid_vel_y(isCoupled("fluid_vel_y") ? coupledValue("fluid_vel_y") : _zero),
    _fluid_vel_y_id(coupled("fluid_vel_y")),
    _fluid_vel_z(isCoupled("fluid_vel_z") ? coupledValue("fluid_vel_z") : _zero),
    _fluid_vel_z_id(coupled("fluid_vel_z")),
    _rho_f(getMaterialPropertyByName<Real>("rho_f")),
    _biot_p_wave_speed(getParam<Real>("p_wave_speed"))
{
}

Real
BiotSlowWaveNonReflectDashpotBC_x::computeQpResidual()
{
  std::vector<Real> velocity(3, 0.0);
  velocity[0] = _u[_qp];
  velocity[1] = _fluid_vel_y[_qp];
  velocity[2] = _fluid_vel_z[_qp];

  Real normal_vel = 0.0;
  for (unsigned int i = 0; i < 3; ++i){
    normal_vel += velocity[i] * _normals[_qp](i);
  }

  return _test[_i][_qp] * _rho_f[_qp] *
         (_biot_p_wave_speed * normal_vel * _normals[_qp](_component));
}
Real
BiotSlowWaveNonReflectDashpotBC_x::computeQpJacobian()
{

  return _test[_i][_qp] * _rho_f[_qp] * _phi[_j][_qp] *
          (_biot_p_wave_speed  * _normals[_qp](_component) * _normals[_qp](_component));
}
Real
BiotSlowWaveNonReflectDashpotBC_x::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _fluid_vel_y_id)
    return _test[_i][_qp] * _rho_f[_qp] *
         (_biot_p_wave_speed * _normals[_qp](_fluid_vel_y_id) * _normals[_qp](_component))* _phi[_j][_qp];
  else if (jvar == _fluid_vel_z_id)
   return _test[_i][_qp] * _rho_f[_qp] *
         (_biot_p_wave_speed * _normals[_qp](_fluid_vel_z_id) * _normals[_qp](_component))* _phi[_j][_qp];
  return 0.0;
}

