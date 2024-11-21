//Implement Lysmer damper ( same as NonreflectingBC | MASTODON ) but in explicit scheme ( follow DashpotBC | MOOSE )
//Created by Amr Ibrahim, April 18, 2024

#include "GeneralNonReflectDashpotBC.h"

registerMooseObject("farmsApp", GeneralNonReflectDashpotBC);

InputParameters
GeneralNonReflectDashpotBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredParam<unsigned int>(
      "component", "The displacement component corresponding the variable this BC acts on.");
  params.addRequiredCoupledVar("displacements", "Displacement in the x,y,z direction");
  params.addRequiredRangeCheckedParam<Real>(
      "p_wave_speed", "p_wave_speed>0.0", "P-wave speed of the material.");
  params.addRequiredRangeCheckedParam<Real>(
      "shear_wave_speed", "shear_wave_speed>0.0", "shear wave speed of the material.");

  return params;
}

GeneralNonReflectDashpotBC::GeneralNonReflectDashpotBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(coupledIndices("displacements")),
    _disp_dot(coupledVectorDot("displacements")),
    _dvel_dot(coupledVectorDotDu("displacements")),
    _du_dot_du(_var.duDotDu()),
    _density(getMaterialPropertyByName<Real>("density")),
    _p_wave_speed(getParam<Real>("p_wave_speed")),
    _shear_wave_speed(getParam<Real>("shear_wave_speed"))
{
}

Real
GeneralNonReflectDashpotBC::computeQpResidual()
{
  Real normal_vel = 0.0;
  for (unsigned int i = 0; i < _ndisp; ++i){
    normal_vel += _disp_dot[_qp](i) * _normals[_qp](i);
  }

  return _test[_i][_qp] * _density[_qp] *
         (_p_wave_speed * normal_vel * _normals[_qp](_component) +
          _shear_wave_speed * ( _disp_dot[_qp](_component) - normal_vel * _normals[_qp](_component)));
}
Real
GeneralNonReflectDashpotBC::computeQpJacobian()
{
  // see PenaltyInclinedNoDisplacementBC.C  
  return _test[_i][_qp] * _density[_qp] *
         (_p_wave_speed *  _du_dot_du[_qp] * _normals[_qp](_component) * _normals[_qp](_component) +
          _shear_wave_speed * ( _du_dot_du[_qp] -  _du_dot_du[_qp] * _normals[_qp](_component) * _normals[_qp](_component)))* _phi[_j][_qp];
}
Real
GeneralNonReflectDashpotBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  // see PenaltyInclinedNoDisplacementBC.C  
  for (unsigned int coupled_component = 0; coupled_component < _ndisp; ++coupled_component)
    if (jvar == _disp_var[coupled_component])
    {
      return _test[_i][_qp] * _density[_qp] *
         (_p_wave_speed * _dvel_dot[_qp]* _normals[_qp](coupled_component) * _normals[_qp](_component) +
          _shear_wave_speed * (_dvel_dot[_qp] - _dvel_dot[_qp]* _normals[_qp](coupled_component) * _normals[_qp](_component)))* _phi[_j][_qp]; 
    }
  return 0.0;
}

