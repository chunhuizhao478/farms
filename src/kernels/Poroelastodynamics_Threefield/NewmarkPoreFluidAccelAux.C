
#include "NewmarkPoreFluidAccelAux.h"

registerMooseObject("farmsApp", NewmarkPoreFluidAccelAux);

InputParameters
NewmarkPoreFluidAccelAux::validParams()
{
InputParameters params = AuxKernel::validParams();
params.addRequiredCoupledVar("darcyvel","Darcy Velocity");
params.addRequiredParam<Real>("gamma","gamma parameter");
return params;
}
NewmarkPoreFluidAccelAux::NewmarkPoreFluidAccelAux(const InputParameters & parameters) :
AuxKernel(parameters),
_w_old(coupledValueOld("darcyvel")),
_w(coupledValue("darcyvel")),
_u_old(valueOld()),
_gamma(getParam<Real>("gamma"))
{
}
Real
NewmarkPoreFluidAccelAux::computeValue()
{
if (!isNodal())
mooseError("NewmarkPoreFluidAccelAux must run on a nodal variable");
Real af_old = _u_old[_qp];
if (_dt == 0)
return af_old;
return 1.0/_gamma*((_w[_qp]-_w_old[_qp])/_dt - af_old*(1.0-_gamma));
}