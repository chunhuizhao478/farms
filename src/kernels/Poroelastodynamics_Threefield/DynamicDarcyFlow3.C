#include "DynamicDarcyFlow3.h"
#include "SubProblem.h"
#include "TimeIntegrator.h"

registerMooseObject("farmsApp", DynamicDarcyFlow3);

InputParameters
DynamicDarcyFlow3::validParams()
{
InputParameters params = TimeKernel::validParams();
params.set<bool>("use_displaced_mesh") = false;
params.addCoupledVar("skeleton_acceleration","skeleton acceleration variable for explicit solver");
return params;
}
DynamicDarcyFlow3::DynamicDarcyFlow3(const InputParameters & parameters)
:TimeKernel(parameters),
_rhof(getMaterialProperty<Real>("rhof")),
_taut(getMaterialProperty<Real>("taut")),
_nf(getMaterialProperty<Real>("porosity")),
_K(getMaterialProperty<Real>("hydconductivity")),
_u_older(valueOlder()),
_us(coupledValue("skeleton_acceleration")),
_us_old(coupledValueOld("skeleton_acceleration")),
_us_older(coupledValueOlder("skeleton_acceleration")),
_us_var_num(coupled("skeleton_acceleration"))

{

}
Real
DynamicDarcyFlow3::computeQpResidual()
{
if (_dt == 0)
    return 0.0;

Real rhoa = (_taut[_qp] - 1) * _nf[_qp];
  
// Central difference for acceleration (second time derivative)
Real us_acc = (_us[_qp] - 2.0 * _us_old[_qp] + _us_older[_qp]) / (_dt * _dt);
  
// Central difference for velocity (first time derivative)
Real u_vel = (_u[_qp] - _u_older[_qp]) / (2.0 * _dt);
  
return _test[_i][_qp] * _rhof[_qp] * (us_acc + u_vel * (1 + rhoa/_nf[_qp])/_nf[_qp]) + 
       _test[_i][_qp] * (1.0/_K[_qp]) * _u[_qp];

}
Real
DynamicDarcyFlow3::computeQpJacobian()
{
if (_dt == 0)
    return 0.0;

Real rhoa = (_taut[_qp] - 1) * _nf[_qp];
  
// Derivative of velocity term with respect to u
Real velocity_jac = 1.0 / (2.0 * _dt) * (1 + rhoa/_nf[_qp])/_nf[_qp];
  
return _test[_i][_qp] * _rhof[_qp] * velocity_jac * _phi[_j][_qp] + 
       _test[_i][_qp] * (1.0/_K[_qp]) * _phi[_j][_qp];

}
Real
DynamicDarcyFlow3::computeQpOffDiagJacobian(unsigned int jvar)
{
if (_dt == 0)
    return 0.0;

if (jvar != _us_var_num)
    return 0.0;

// Derivative of acceleration term with respect to us
Real acc_jac = 1.0 / (_dt * _dt);
  
return _test[_i][_qp] * _rhof[_qp] * acc_jac * _phi[_j][_qp];

}
