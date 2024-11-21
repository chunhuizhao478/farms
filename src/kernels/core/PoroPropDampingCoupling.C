/*
Kernel of Stiffness Proportional Damping 
*/

#include "PoroPropDampingCoupling.h"

registerMooseObject("farmsApp", PoroPropDampingCoupling);

InputParameters
PoroPropDampingCoupling::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription(
      "Compute Stiffness Proportional Damping Residual");
  params.addRequiredCoupledVar("porepressure", "Pore pressure, $p_s$.");
  params.addRequiredParam<unsigned int>("component",
                                        "The gradient direction (0 for x, 1 for y and 2 for z)");
  params.addRequiredParam<Real>("q","Ratio Factor to assign magnitude of stiffness proportional damping term");

  return params;
}

PoroPropDampingCoupling::PoroPropDampingCoupling(const InputParameters & parameters)
  :Kernel(parameters),
  _component(getParam<unsigned int>("component")),
  _coefficient(getMaterialProperty<Real>("biot_coefficient")),
  _porepressure(coupledValue("porepressure")),
  _porepressure_old(coupledValueOlder("porepressure")),
  _porepressure_var_num(coupled("porepressure")),
  //Ratio factor 
  _q(getParam<Real>("q"))

{
}

Real
PoroPropDampingCoupling::computeQpResidual()
{
  return -_coefficient[_qp] * ( _porepressure[_qp] - _porepressure_old[_qp] ) * _grad_test[_i][_qp](_component) * _q;
}

Real
PoroPropDampingCoupling::computeQpJacobian()
{
  if (_var.number() != _porepressure_var_num)
  {
    return 0.0;
  }
  else
  {
    return -_coefficient[_qp] * _phi[_j][_qp] * _grad_test[_i][_qp](_component) * _q;
  }
}

Real
PoroPropDampingCoupling::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar != _porepressure_var_num)
  {
    return 0.0;
  }
  else
  {
    return -_coefficient[_qp] * _phi[_j][_qp] * _grad_test[_i][_qp](_component) *_q;
  }
    
  
}



