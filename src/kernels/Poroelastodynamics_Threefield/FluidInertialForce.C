#include "FluidInertialForce.h"
#include "SubProblem.h"
#include "TimeIntegrator.h"

registerMooseObject("farmsApp", FluidInertialForce);

InputParameters
FluidInertialForce::validParams()
{
InputParameters params = TimeKernel::validParams();
params.set<bool>("use_displaced_mesh") = false;
params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");
return params;
}
FluidInertialForce::FluidInertialForce(const InputParameters & parameters)
:TimeKernel(parameters),
_lumping(getParam<bool>("lumping")),
_rhof(getMaterialProperty<Real>("rhof")),
_nf(getMaterialProperty<Real>("porosity")),
_time_integrator(*_sys.getTimeIntegrator())
{
    _du_dot_du = &(_var.duDotDu());


  if (_time_integrator.isLumped())
    {
    _u_dot_factor_dof = &_var.vectorTagDofValue(_time_integrator.uDotFactorTag());
    }
    _u_dot_factor = &_var.vectorTagValue(_time_integrator.uDotFactorTag());
}
Real
FluidInertialForce::computeQpResidual()
{
if (_dt == 0)
return 0.0;
return _test[_i][_qp]*_rhof[_qp]* (*_u_dot_factor)[_qp]/_nf[_qp];
if (_time_integrator.isLumped())
  {
  return _test[_i][_qp] * _rhof[_qp] /_nf[_qp] ;
  for (unsigned int i = 0; i < _test.size(); ++i)
  this->_local_re(i) *=  (*_u_dot_factor_dof)[i];
  }
}
Real
FluidInertialForce::computeQpJacobian()
{
if (_dt == 0)
return 0.0;
return _test[_i][_qp]*_rhof[_qp]* (*_du_dot_du)[_qp] /_nf[_qp] ;
}