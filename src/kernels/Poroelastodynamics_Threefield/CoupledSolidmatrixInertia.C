#include "CoupledSolidmatrixInertia.h"
#include "SubProblem.h"
#include "TimeIntegrator.h"

registerMooseObject("farmsApp", CoupledSolidmatrixInertia);

InputParameters
CoupledSolidmatrixInertia::validParams()
{
InputParameters params = TimeKernel::validParams();
params.set<bool>("use_displaced_mesh") = false;
params.addCoupledVar("displacements","skeleton acceleration variable for explicit solver");
params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");
return params;
}
CoupledSolidmatrixInertia::CoupledSolidmatrixInertia(const InputParameters & parameters)
:TimeKernel(parameters),
_lumping(getParam<bool>("lumping")),
_rhof(getMaterialProperty<Real>("rhof")),
_us_var_num(coupled("displacements")),
_us_older(coupledDofValuesOlder("displacements")),
_us_old(coupledDofValuesOld("displacements")),
_us(coupledDofValues("displacements")),
_time_integrator(_sys.getTimeIntegrator())
{
    this->addFEVariableCoupleableVectorTag(_time_integrator.uDotFactorTag());
    this->addFEVariableCoupleableVectorTag(_time_integrator.uDotDotFactorTag());

    _us_dot_dot = &coupledVectorTagValue("displacements",_time_integrator.uDotDotFactorTag());
    _us_dot_dot_dof = &coupledVectorTagDofValue("displacements",_time_integrator.uDotDotFactorTag());

    _dus_dot_dot_du= &this->coupledDotDotDu("displacements");

     _u_dot_factor = &_var.vectorTagValue(_time_integrator.uDotFactorTag());
      _du_dot_du = &(_var.duDotDu());

}
Real
CoupledSolidmatrixInertia::computeQpResidual()
{
if (_dt == 0)
return 0.0;
return _test[_i][_qp]*_rhof[_qp]*( _us[_i] - 2 *_us_old[_i] + _us_older[_i])/_dt/_dt;
 if (_time_integrator.isLumped())
 {
 return _test[_i][_qp] * _rhof[_qp] /_dt/_dt;
 for (unsigned int i = 0; i < _test.size(); ++i)
 this->_local_re(i) *=  ( _us[_i] - 2 *_us_old[_i] + _us_older[_i]);
  }
}
Real
CoupledSolidmatrixInertia::computeQpJacobian()
{
return 0.0;
}
//return _test[_i][_qp]*_rhof[_qp]* ((*_du_dot_du)[_qp] / _nf[_qp]+_gravity/_K[_qp])*_phi[_j][_qp];

Real
CoupledSolidmatrixInertia::computeQpOffDiagJacobian(unsigned int jvar)
{
if (_dt == 0)
return 0.0;
if (jvar != _us_var_num) 
return 0.0;
return _test[_i][_qp]*_rhof[_qp]*(*_dus_dot_dot_du)[_qp]*_phi[_j][_qp];
}