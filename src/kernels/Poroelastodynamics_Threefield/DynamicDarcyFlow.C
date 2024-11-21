#include "DynamicDarcyFlow.h"
#include "SubProblem.h"

registerMooseObject("farmsApp", DynamicDarcyFlow);

InputParameters
DynamicDarcyFlow::validParams()
{
InputParameters params = Kernel::validParams();
params.set<bool>("use_displaced_mesh") = false;
params.addRequiredCoupledVar("skeletondisp","skeleton displacement variable");
params.addRequiredCoupledVar("skeletonvel","skeleton velocity variable");
params.addRequiredCoupledVar("skeletonaccel","skeleton acceleration variable");
params.addRequiredCoupledVar("fluidaccel","fluid relative acceleration variable");
params.addRequiredParam<Real>("gravity","acceleration due to gravity");
params.addRequiredParam<Real>("beta","beta parameter");
params.addRequiredParam<Real>("gamma","gamma parameter");
return params;
}
DynamicDarcyFlow::DynamicDarcyFlow(const InputParameters & parameters)
:Kernel(parameters),
_rhof(getMaterialProperty<Real>("rhof")),
_nf(getMaterialProperty<Real>("porosity")),
_K(getMaterialProperty<Real>("hydconductivity")),
_us(coupledValue("skeletondisp")),
_us_old(coupledValueOld("skeletondisp")),
_vs_old(coupledValueOld("skeletonvel")),
_as_old(coupledValueOld("skeletonaccel")),
_u_old(valueOld()),
_af_old(coupledValueOld("fluidaccel")),
_us_var_num(coupled("skeletondisp")),
_gravity(getParam<Real>("gravity")),
_beta(getParam<Real>("beta")),
_gamma(getParam<Real>("gamma"))
{}
Real
DynamicDarcyFlow::computeQpResidual()
{
if (_dt == 0)
return 0.0;
Real as=1/_beta*(((_us[_qp]-_us_old[_qp])/(_dt*_dt)) - _vs_old[_qp]/_dt -
_as_old[_qp]*(0.5-_beta));
Real af=1/_gamma*((_u[_qp]-_u_old[_qp])/_dt - (1.0-_gamma)*_af_old[_qp]);
return _test[_i][_qp]*_rhof[_qp]*( as + af/_nf[_qp] + _gravity/_K[_qp]*_u[_qp] );
}
Real
DynamicDarcyFlow::computeQpJacobian()
{
if (_dt == 0)
return 0.0;
return _test[_i][_qp]*_rhof[_qp]*(1.0/(_nf[_qp]*_gamma*_dt) +
_gravity/_K[_qp])*_phi[_j][_qp];
}
Real
DynamicDarcyFlow::computeQpOffDiagJacobian(unsigned int jvar)
{
if (_dt == 0)
return 0.0;
if (jvar != _us_var_num) // only u-w block is nonzero
return 0.0;
return _test[_i][_qp]*_rhof[_qp]/(_beta*_dt*_dt)*_phi[_j][_qp];
}