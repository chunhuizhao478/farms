#include "INSmassSolid.h"
#include "TimeIntegrator.h"

registerMooseObject("farmsApp", INSmassSolid);

InputParameters
INSmassSolid::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.set<bool>("use_displaced_mesh") = false;
  params.addCoupledVar("displacements", "String of displacement components");
  return params;
}

INSmassSolid::INSmassSolid(const InputParameters & parameters)
  : TimeKernel(parameters),
    _coefficient(getMaterialProperty<Real>("biot_coefficient")),
    _ux_var(coupled("displacements", 0)),
    _uy_var(_mesh.dimension() >= 2 ? coupled("displacements", 1) : libMesh::invalid_uint),
    _uz_var(_mesh.dimension() == 3 ? coupled("displacements", 2) : libMesh::invalid_uint),
    _vx_var(&getVar("displacements", 0)),
    _vy_var(&getVar("displacements", 1)),
    _vz_var(&getVar("displacements", 2)),
    _vectorx_phi(_assembly.gradPhi(*_vx_var)),
    _vectory_phi(_assembly.gradPhi(*_vy_var)),
    _vectorz_phi(_assembly.gradPhi(*_vz_var)),
    _time_integrator(_sys.getTimeIntegrator(_var.number()))  // Fixed: pass variable number
{
  this->addFEVariableCoupleableVectorTag(_time_integrator.uDotFactorTag());
  _grad_ux_dot = &coupledVectorTagGradient("displacements", _time_integrator.uDotFactorTag(), 0);
  _grad_uy_dot = &coupledVectorTagGradient("displacements", _time_integrator.uDotFactorTag(), 1);
  
  if (_mesh.dimension() == 3)
    _grad_uz_dot = &coupledVectorTagGradient("displacements", _time_integrator.uDotFactorTag(), 2);
    
  _d_grad_ux_dot_dv = &this->coupledDotDu("displacements", 0);
  _d_grad_uy_dot_dv = &this->coupledDotDu("displacements", 1);
  
  if (_mesh.dimension() == 3)
    _d_grad_uz_dot_dv = &this->coupledDotDu("displacements", 2);
}

Real
INSmassSolid::computeQpResidual()
{
  Real div_u_dot = (*_grad_ux_dot)[_qp](0) + (*_grad_uy_dot)[_qp](1);
  
  if (_mesh.dimension() == 3)
    div_u_dot += (*_grad_uz_dot)[_qp](2);
    
  return _coefficient[_qp] * _test[_i][_qp] * div_u_dot;
}

Real
INSmassSolid::computeQpJacobian()
{
  return 0.0;
}

Real
INSmassSolid::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_dt == 0)
    return 0.0;
    
  if (jvar == _ux_var)
    return _vectorx_phi[_j][_qp](0) * _test[_i][_qp] * _coefficient[_qp] * (*_d_grad_ux_dot_dv)[_qp];
  else if (jvar == _uy_var)
    return _vectory_phi[_j][_qp](1) * _test[_i][_qp] * _coefficient[_qp] * (*_d_grad_uy_dot_dv)[_qp];
  else if (_mesh.dimension() == 3 && jvar == _uz_var)
    return _vectorz_phi[_j][_qp](2) * _test[_i][_qp] * _coefficient[_qp] * (*_d_grad_uz_dot_dv)[_qp];
  else
    return 0.0;
}