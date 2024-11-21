#include "DarcyFlow.h"
#include "SubProblem.h"
#include "TimeIntegrator.h"

registerMooseObject("farmsApp", DarcyFlow);

InputParameters
DarcyFlow::validParams()
{
InputParameters params = TimeKernel::validParams();
params.set<bool>("use_displaced_mesh") = false;
return params;
}
DarcyFlow::DarcyFlow(const InputParameters & parameters)
:TimeKernel(parameters),
_K(getMaterialProperty<Real>("hydconductivity"))
{

}
Real
DarcyFlow::computeQpResidual()
{
if (_dt == 0)
return 0.0;
return _test[_i][_qp]*1/_K[_qp]*_u[_qp] ;
}
Real
DarcyFlow::computeQpJacobian()
{
if (_dt == 0)
return 0.0;
return _test[_i][_qp]*1/_K[_qp]*_phi[_j][_qp];
}
