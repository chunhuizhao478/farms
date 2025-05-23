#include "FluidStorageNewmark.h"

registerMooseObject("farmsApp", FluidStorageNewmark);

InputParameters
FluidStorageNewmark::validParams()
{
InputParameters params = Kernel::validParams();
params.set<bool>("use_displaced_mesh") = false;
return params;
}
FluidStorageNewmark::FluidStorageNewmark(const InputParameters & parameters)
:Kernel(parameters),
_u_old(valueOld()),
_coefficient(getMaterialProperty<Real>("biot_modulus"))
{}
Real
FluidStorageNewmark::computeQpResidual()

{
if (_dt == 0)
return 0.0;

Real p_dot = (_u[_qp]-_u_old[_qp])/_dt;

return _test[_i][_qp] * p_dot *  1/_coefficient[_qp];
}
Real
FluidStorageNewmark::computeQpJacobian()
{
// Derivative wrt p 
if (_dt == 0)
return 0.0;
return  _test[_i][_qp] *_phi[_j][_qp] /_dt * 1/_coefficient[_qp];
}
Real
FluidStorageNewmark::computeQpOffDiagJacobian(unsigned int jvar)
{
if (_dt == 0)
return 0.0;

return 0.0;
}