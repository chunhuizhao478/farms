#include "DynamicDarcyFlow2.h"
#include "SubProblem.h"
#include "TimeIntegrator.h"

registerMooseObject("farmsApp", DynamicDarcyFlow2);

InputParameters
DynamicDarcyFlow2::validParams()
{
InputParameters params = TimeKernel::validParams();
params.set<bool>("use_displaced_mesh") = false;
params.addCoupledVar("skeleton_acceleration","skeleton acceleration variable for explicit solver");
params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");
return params;
}
DynamicDarcyFlow2::DynamicDarcyFlow2(const InputParameters & parameters)
:TimeKernel(parameters),
_lumping(getParam<bool>("lumping")),
_rhof(getMaterialProperty<Real>("rhof")),
_taut(getMaterialProperty<Real>("taut")),
_nf(getMaterialProperty<Real>("porosity")),
_K(getMaterialProperty<Real>("hydconductivity")),
//_us_dot_dot(coupledDotDot("skeleton_acceleration")),
_us_var_num(coupled("skeleton_acceleration")),
_time_integrator(*_sys.getTimeIntegrator())
{
//    _u_dot_old = &(_var.uDotOld());
    _du_dot_du = &(_var.duDotDu());

    this->addFEVariableCoupleableVectorTag(_time_integrator.uDotFactorTag());
    this->addFEVariableCoupleableVectorTag(_time_integrator.uDotDotFactorTag());

    _us_dot_dot = &coupledVectorTagValue("skeleton_acceleration",_time_integrator.uDotDotFactorTag());

    _dus_dot_dot_du= &this->coupledDotDotDu("skeleton_acceleration");

  if (_time_integrator.isLumped())
    {
    _u_dot_factor_dof = &_var.vectorTagDofValue(_time_integrator.uDotFactorTag());
    }
    _u_dot_factor = &_var.vectorTagValue(_time_integrator.uDotFactorTag());
}
Real
DynamicDarcyFlow2::computeQpResidual()
{
if (_dt == 0)
return 0.0;

Real rhoa = (_taut[_qp]-1)*_nf[_qp];

return _test[_i][_qp]*_rhof[_qp]*( (*_us_dot_dot)[_qp] + (*_u_dot_factor)[_qp]*(1+rhoa/_nf[_qp])/_nf[_qp] + 1/_K[_qp]/_rhof[_qp]*_u[_qp] );


}
Real
DynamicDarcyFlow2::computeQpJacobian()
{
if (_dt == 0)
return 0.0;
Real rhoa = (_taut[_qp]-1)*_nf[_qp];
return _test[_i][_qp]*_rhof[_qp]*( (*_du_dot_du)[_qp] *(1+rhoa/_nf[_qp])/_nf[_qp] +
1/_K[_qp]/_rhof[_qp])*_phi[_j][_qp];

}
Real
DynamicDarcyFlow2::computeQpOffDiagJacobian(unsigned int jvar)
{
if (_dt == 0)
return 0.0;
if (jvar != _us_var_num) // only u-w block is nonzero
return 0.0;
return _test[_i][_qp]*_rhof[_qp]*(*_dus_dot_dot_du)[_qp]*_phi[_j][_qp];

}
void
DynamicDarcyFlow2::computeJacobian()
{
  if (_lumping)
  {
    prepareMatrixTag(_assembly, _var.number(), _var.number());

    precalculateJacobian();
    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < _phi.size(); _j++)
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
          _local_ke(_i, _i) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();

    accumulateTaggedLocalMatrix();
  }
  else
    TimeKernel::computeJacobian();
}