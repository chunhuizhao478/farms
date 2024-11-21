#include "CoupledFluidInertialForce.h"
#include "SubProblem.h"
#include "TimeIntegrator.h"

registerMooseObject("farmsApp", CoupledFluidInertialForce);

InputParameters
CoupledFluidInertialForce::validParams()
{
InputParameters params = TimeKernel::validParams();
params.set<bool>("use_displaced_mesh") = false;
params.addCoupledVar("fluid_vel","fluid relative acceleration variable for explicit solver");
params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");
return params;
}
CoupledFluidInertialForce::CoupledFluidInertialForce(const InputParameters & parameters)
:TimeKernel(parameters),
_lumping(getParam<bool>("lumping")),
_rhof(getMaterialProperty<Real>("rhof")),
_w_var_num(coupled("fluid_vel")),
_wf(coupledDofValues("fluid_vel")),
_wf_older(coupledDofValuesOlder("fluid_vel")),
_time_integrator(*_sys.getTimeIntegrator())
{

    this->addFEVariableCoupleableVectorTag(_time_integrator.uDotFactorTag());

    _wf_dot_factor_dof = &coupledVectorTagDofValue("fluid_vel",_time_integrator.uDotFactorTag());
    _wf_dot_factor = &coupledVectorTagValue("fluid_vel",_time_integrator.uDotFactorTag());

    _dwf_dot_du= &this->coupledDotDu("fluid_vel");
}
Real
CoupledFluidInertialForce::computeQpResidual()
{
if (_dt == 0)
return 0.0;
return _test[_i][_qp] * _rhof[_qp] * (*_wf_dot_factor)[_qp];
 if (_time_integrator.isLumped())
  {
   return _test[_i][_qp] * _rhof[_qp] /2/_dt  ;
  for (unsigned int i = 0; i < _test.size(); ++i)
   this->_local_re(i) *= (_wf[_i]-_wf_older[_i]);
   }
}
Real
CoupledFluidInertialForce::computeQpJacobian()
{
return 0.0;
}
Real
CoupledFluidInertialForce::computeQpOffDiagJacobian(unsigned int jvar)
{
if (_dt == 0)
return 0.0;
if (jvar != _w_var_num) 
return 0.0;
return _test[_i][_qp]*_rhof[_qp]*(*_dwf_dot_du)[_qp]*_phi[_j][_qp];
}