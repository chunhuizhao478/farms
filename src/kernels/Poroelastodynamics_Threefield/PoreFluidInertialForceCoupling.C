#include "PoreFluidInertialForceCoupling.h"
#include "SubProblem.h"

registerMooseObject("farmsApp", PoreFluidInertialForceCoupling);

InputParameters
PoreFluidInertialForceCoupling::validParams()
{
InputParameters params = Kernel::validParams();
params.set<bool>("use_displaced_mesh") = false;
params.addRequiredCoupledVar("fluidaccel","fluid relative acceleration variable");
params.addRequiredCoupledVar("darcyvel","Darcy velocity variable");
params.addRequiredParam<Real>("gamma","gamma parameter");
return params;
}
PoreFluidInertialForceCoupling::PoreFluidInertialForceCoupling(const InputParameters & parameters)
:Kernel(parameters),
_rhof(getMaterialProperty<Real>("rhof")),
_af_old(coupledValueOld("fluidaccel")),
_w(coupledValue("darcyvel")),
_w_old(coupledValueOld("darcyvel")),
_w_var_num(coupled("darcyvel")),
_gamma(getParam<Real>("gamma"))
{}
Real
PoreFluidInertialForceCoupling::computeQpResidual()
{
if (_dt == 0)
return 0.0;
Real af=1/_gamma*((_w[_qp] - _w_old[_qp])/_dt - (1.0-_gamma)*_af_old[_qp]);
return _test[_i][_qp]*_rhof[_qp]*af;
}
Real
PoreFluidInertialForceCoupling::computeQpJacobian()
{
return 0.0;
}
Real
PoreFluidInertialForceCoupling::computeQpOffDiagJacobian(unsigned int jvar)
{
if (_dt == 0)
return 0.0;
if (jvar != _w_var_num) // only u-w block is nonzero
return 0.0;
return _test[_i][_qp]*_rhof[_qp]/(_gamma*_dt)*_phi[_j][_qp];
}