#include "FluidStorageNewmark.h"

registerMooseObject("farmsApp", FluidStorageNewmark);

InputParameters
FluidStorageNewmark::validParams()
{
InputParameters params = Kernel::validParams();
params.set<bool>("use_displaced_mesh") = false;
params.addRequiredParam<Real>("gamma","gamma parameter");
params.addRequiredCoupledVar("pressure_dot","time derivative of pressure variable");
return params;
}
FluidStorageNewmark::FluidStorageNewmark(const InputParameters & parameters)
:Kernel(parameters),
_gamma(getParam<Real>("gamma")),
_u_old(valueOld()),
_p_dot_old(coupledValueOld("pressure_dot")),
_coefficient(getMaterialProperty<Real>("biot_modulus"))
{}
Real
FluidStorageNewmark::computeQpResidual()

{
if (_dt == 0)
return 0.0;

Real p_dot = 1.0/_gamma*((_u[_qp]-_u_old[_qp])/_dt - _p_dot_old[_qp]*(1.0-_gamma));

return -_test[_i][_qp] * p_dot *  1/_coefficient[_qp];
}
Real
FluidStorageNewmark::computeQpJacobian()
{
// Derivative wrt p 
if (_dt == 0)
return 0.0;
return - _test[_i][_qp] * 1.0/(_gamma*_dt) *_phi[_j][_qp] * 1/_coefficient[_qp];
}
Real
FluidStorageNewmark::computeQpOffDiagJacobian(unsigned int jvar)
{
if (_dt == 0)
return 0.0;

return 0.0;
}