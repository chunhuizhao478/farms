//Implement Lysmer damper ( same as NonreflectingBC | MASTODON ) but in explicit scheme ( follow DashpotBC | MOOSE )
//Generationalization to 3D
//Created by Chunhui Zhao, May 1, 2023

#include "NonReflectDashpotBC3d.h"

registerMooseObject("farmsApp", NonReflectDashpotBC3d);

InputParameters
NonReflectDashpotBC3d::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredParam<unsigned int>(
      "component", "The displacement component corresponding the variable this BC acts on.");
  params.addRequiredCoupledVar("disp_x", "Displacement in the x direction");
  params.addRequiredCoupledVar("disp_y", "Displacement in the y direction");
  params.addRequiredCoupledVar("disp_z", "Displacement in the z direction");
  params.addRequiredRangeCheckedParam<Real>(
      "p_wave_speed", "p_wave_speed>0.0", "P-wave speed of the material.");
  params.addRequiredRangeCheckedParam<Real>(
      "shear_wave_speed", "shear_wave_speed>0.0", "shear wave speed of the material.");

  return params;
}

NonReflectDashpotBC3d::NonReflectDashpotBC3d(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _component(getParam<unsigned int>("component")),
    _disp_x_var(coupled("disp_x")),
    _disp_y_var(coupled("disp_y")),
    _disp_z_var(coupled("disp_z")),
    _disp_x_dot(coupledDot("disp_x")),
    _disp_y_dot(coupledDot("disp_y")),
    _disp_z_dot(coupledDot("disp_z")),
    _density(getMaterialPropertyByName<Real>("density")),
    _p_wave_speed(getParam<Real>("p_wave_speed")),
    _shear_wave_speed(getParam<Real>("shear_wave_speed"))
{
}

Real
NonReflectDashpotBC3d::computeQpResidual()
{
  std::vector<Real> velocity(3, 0.0);
  velocity[0] = _disp_x_dot[_qp];
  velocity[1] = _disp_y_dot[_qp];
  velocity[2] = _disp_z_dot[_qp];

  Real normal_vel = 0.0;
  for (unsigned int i = 0; i < 3; ++i){
    normal_vel += velocity[i] * _normals[_qp](i);
  }

  return _test[_i][_qp] * _density[_qp] *
         (_p_wave_speed * normal_vel * _normals[_qp](_component) +
          _shear_wave_speed * (velocity[_component] - normal_vel * _normals[_qp](_component)));
}