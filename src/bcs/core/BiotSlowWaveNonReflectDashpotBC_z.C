//Implement Lysmer damper ( same as NonreflectingBC | MASTODON ) but in explicit scheme ( follow DashpotBC | MOOSE )
//Created by Amr Ibrahim, April 18, 2024

#include "BiotSlowWaveNonReflectDashpotBC_z.h"

registerMooseObject("farmsApp", BiotSlowWaveNonReflectDashpotBC_z);

InputParameters
BiotSlowWaveNonReflectDashpotBC_z::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredParam<unsigned int>(
      "component", "The Fluid velocity component corresponding the variable this BC acts on.");
  params.addCoupledVar("fluid_vel_x", "Fluid velocity in the y direction");
  params.addCoupledVar("fluid_vel_y", "Fluid velocity in the Z direction");
  params.addRequiredRangeCheckedParam<Real>(
      "_biot_p_wave_speed", "p_wave_speed>0.0", "P-wave speed of the material.");
  return params;
}

BiotSlowWaveNonReflectDashpotBC_z::BiotSlowWaveNonReflectDashpotBC_z(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _fluid_vel_x(isCoupled("fluid_vel_x") ? coupledValue("fluid_vel_x") : _zero),
    _fluid_vel_x_id(coupled("fluid_vel_x")),
    _fluid_vel_y(isCoupled("fluid_vel_y") ? coupledValue("fluid_vel_y") : _zero),
    _fluid_vel_y_id(coupled("fluid_vel_y")),
    _rho_f(getMaterialPropertyByName<Real>("rho_f")),
    _biot_p_wave_speed(getParam<Real>("p_wave_speed"))
{
}

Real
BiotSlowWaveNonReflectDashpotBC_z::computeQpResidual()
{
  std::vector<Real> velocity(3, 0.0);
  velocity[0] = _fluid_vel_x[_qp];
  velocity[1] = _fluid_vel_y[_qp];
  velocity[2] = _u[_qp];

  Real normal_vel = 0.0;
  for (unsigned int i = 0; i < 3; ++i){
    normal_vel += velocity[i] * _normals[_qp](i);
  }

  return _test[_i][_qp] * _rho_f[_qp] *
         (_biot_p_wave_speed * normal_vel * _normals[_qp](_component));
}
Real
BiotSlowWaveNonReflectDashpotBC_z::computeQpJacobian()
{
  // see PenaltyInclinedNoDisplacementBC.C  
  return _test[_i][_qp] * _rho_f[_qp] * _phi[_j][_qp] *
         (_biot_p_wave_speed  * _normals[_qp](_component) * _normals[_qp](_component));
}
Real
BiotSlowWaveNonReflectDashpotBC_z::computeQpOffDiagJacobian(unsigned int jvar)
{
  // see PenaltyInclinedNoDisplacementBC.C  
  if (jvar == _fluid_vel_x_id)
    return _test[_i][_qp] * _rho_f[_qp] *
         (_biot_p_wave_speed * _normals[_qp](_fluid_vel_x_id) * _normals[_qp](_component))* _phi[_j][_qp];
  else if (jvar == _fluid_vel_y_id)
   return _test[_i][_qp] * _rho_f[_qp] *
         (_biot_p_wave_speed * _normals[_qp](_fluid_vel_y_id) * _normals[_qp](_component))* _phi[_j][_qp];
  return 0.0;
}

