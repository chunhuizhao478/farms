#include "MassConservationNewmark.h"

registerMooseObject("farmsApp", MassConservationNewmark);

InputParameters
MassConservationNewmark::validParams()
{
InputParameters params = Kernel::validParams();
params.set<bool>("use_displaced_mesh") = false;
params.addCoupledVar("displacements", "String of displacement components");
params.addCoupledVar("velocities", "String of velocity components");
params.addCoupledVar("accelerations", "String of acceleration components");
params.addRequiredParam<Real>("beta","beta parameter");
params.addRequiredParam<Real>("gamma","gamma parameter");
return params;
}
MassConservationNewmark::MassConservationNewmark(const InputParameters & parameters)
:Kernel(parameters),
_ux_var(coupled("displacements",0)),
_uy_var(_mesh.dimension() >= 2 ? coupled("displacements",1) : libMesh::invalid_uint),
_uz_var(_mesh.dimension() == 3 ? coupled("displacements",2) : libMesh::invalid_uint),
_grad_ux(coupledGradient("displacements",0)),
_grad_uy(_mesh.dimension() >= 2 ? coupledGradient("displacements",1) : _grad_zero),
_grad_uz(_mesh.dimension() == 3 ? coupledGradient("displacements",2) : _grad_zero),
_grad_ux_old(coupledGradientOld("displacements",0)),
_grad_uy_old(_mesh.dimension() >= 2 ? coupledGradientOld("displacements",1) :
_grad_zero),
_grad_uz_old(_mesh.dimension() == 3 ? coupledGradientOld("displacements",2) :
_grad_zero),
_grad_vx_old(coupledGradientOld("velocities",0)),
_grad_vy_old(_mesh.dimension() >= 2 ? coupledGradientOld("velocities",1) : _grad_zero),
_grad_vz_old(_mesh.dimension() == 3 ? coupledGradientOld("velocities",2) : _grad_zero),
_grad_ax_old(coupledGradientOld("accelerations",0)),
_grad_ay_old(_mesh.dimension() >= 2 ? coupledGradientOld("accelerations",1) :
_grad_zero),
_grad_az_old(_mesh.dimension() == 3 ? coupledGradientOld("accelerations",2) :
_grad_zero),
_beta(getParam<Real>("beta")),
_gamma(getParam<Real>("gamma")),
_coefficient(getMaterialProperty<Real>("biot_coefficient"))
{}
Real
MassConservationNewmark::computeQpResidual()

{
if (_dt == 0)
return 0.0;
Real div_u = _grad_ux[_qp](0) + _grad_uy[_qp](1) + _grad_uz[_qp](2);
Real div_u_old = _grad_ux_old[_qp](0) + _grad_uy_old[_qp](1) + _grad_uz_old[_qp](2);
Real div_v_old = _grad_vx_old[_qp](0) + _grad_vy_old[_qp](1) + _grad_vz_old[_qp](2);
Real div_a_old = _grad_ax_old[_qp](0) + _grad_ay_old[_qp](1) + _grad_az_old[_qp](2);
Real div_v = (_gamma/_beta/_dt)*(div_u - div_u_old) + (1.0-_gamma/_beta)*div_v_old
- (1.0 - _gamma/2.0/_beta)*div_a_old;
return -_test[_i][_qp]*div_v * _coefficient[_qp];
}
Real
MassConservationNewmark::computeQpJacobian()
{
// Derivative wrt p is zero
return 0.0;
}
Real
MassConservationNewmark::computeQpOffDiagJacobian(unsigned int jvar)
{
if (_dt == 0)
return 0.0;
else if (jvar == _ux_var)
return -(_gamma/_beta/_dt)*_grad_phi[_j][_qp](0)*_test[_i][_qp]*_coefficient[_qp];
else if (jvar == _uy_var)
return -(_gamma/_beta/_dt)*_grad_phi[_j][_qp](1)*_test[_i][_qp]*_coefficient[_qp];
else if (jvar == _uz_var)
return -(_gamma/_beta/_dt)*_grad_phi[_j][_qp](2)*_test[_i][_qp]*_coefficient[_qp];
else
return 0.0;
}