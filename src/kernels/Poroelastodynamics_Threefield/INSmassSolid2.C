#include "INSmassSolid2.h"

registerMooseObject("farmsApp", INSmassSolid2);

InputParameters
INSmassSolid2::validParams()
{
InputParameters params = Kernel::validParams();
params.set<bool>("use_displaced_mesh") = false;
params.addCoupledVar("displacements", "String of displacement components");
return params;
}
INSmassSolid2::INSmassSolid2(const InputParameters & parameters)
:Kernel(parameters),
_u_dot(_var.uDot()), 
_du_dot_du(_var.duDotDu()),
_coefficient_M(getMaterialProperty<Real>("biot_modulus")),
_ux_var(coupled("displacements",0)),
_uy_var(_mesh.dimension() >= 2 ? coupled("displacements",1) : libMesh::invalid_uint),
_uz_var(_mesh.dimension() == 3 ? coupled("displacements",2) : libMesh::invalid_uint),
_grad_ux(coupledGradient("displacements",0)),
_grad_uy(_mesh.dimension() >= 2 ? coupledGradient("displacements",1) : _grad_zero),
_grad_uz(_mesh.dimension() == 3 ? coupledGradient("displacements",2) : _grad_zero),
_grad_ux_older(coupledGradientOld("displacements",0)),
_grad_uy_older(_mesh.dimension() >= 2 ? coupledGradientOld("displacements",1) :_grad_zero),
_grad_uz_older(_mesh.dimension() == 3 ? coupledGradientOld("displacements",2) :_grad_zero),
_coefficient(getMaterialProperty<Real>("biot_coefficient"))
{}
Real
INSmassSolid2::computeQpResidual()
{
if (_dt == 0)
return 0.0;

Real div_u = _grad_ux[_qp](0) + _grad_uy[_qp](1) + _grad_uz[_qp](2);
Real div_u_older = _grad_ux_older[_qp](0) + _grad_uy_older[_qp](1) + _grad_uz_older[_qp](2);
Real div_v = ( div_u - div_u_older ) / _dt ;

return _test[_i][_qp] * (_coefficient[_qp] * div_v + _u_dot[_qp] * 1/_coefficient_M[_qp]);
}
Real
INSmassSolid2::computeQpJacobian()
{
// Derivative wrt p is zero
return _test[_i][_qp] * _phi[_j][_qp] * _du_dot_du[_qp] * 1/_coefficient_M[_qp];
}
Real
INSmassSolid2::computeQpOffDiagJacobian(unsigned int jvar)
{
if (_dt == 0)
return 0.0;
else if (jvar == _ux_var)
return ( 1 /_dt) *_grad_phi[_j][_qp](0)*_test[_i][_qp]*_coefficient[_qp];
else if (jvar == _uy_var)
return ( 1 /_dt) *_grad_phi[_j][_qp](1)*_test[_i][_qp]*_coefficient[_qp];
else if (jvar == _uz_var)
return ( 1 /_dt) *_grad_phi[_j][_qp](2)*_test[_i][_qp]*_coefficient[_qp];
else
return 0.0;
}