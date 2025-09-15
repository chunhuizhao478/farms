/*
Kernel of Stiffness Proportional Damping 
*/

#include "LagrangianStiffPropDampingImplicit.h"
#include "ElasticityTensorTools.h"

registerMooseObject("farmsApp", LagrangianStiffPropDampingImplicit);

InputParameters
LagrangianStiffPropDampingImplicit::validParams()
{
  InputParameters params = TotalLagrangianStressDivergence::validParams();
  params.addClassDescription(
      "Compute Stiffness Proportional Damping Residual");
  params.addRequiredRangeCheckedParam<unsigned int>("component",
                                                    "component < 3",
                                                    "An integer corresponding to the direction "
                                                    "the variable this kernel acts in. (0 for x, "
                                                    "1 for y, 2 for z)");
  params.addRequiredParam<Real>("zeta","Ratio Factor to assign magnitude of stiffness proportional damping term");

  return params;
}

LagrangianStiffPropDampingImplicit::LagrangianStiffPropDampingImplicit(const InputParameters & parameters)
  : TotalLagrangianStressDivergence(parameters),
  
  //Get stress tensor from previous time step
  _stress_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "pk1_stress")),
  
  //Get stress tensor from current time step
  _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "pk1_stress")),
  
  //Ratio factor 
  _q(getParam<Real>("zeta")),

  //component
  _component(getParam<unsigned int>("component"))

{
}

Real
LagrangianStiffPropDampingImplicit::computeQpResidual()
{
  if ( _t < _dt)
    return 0.0;
  else
    return _q*(_stress[_qp].row(_component)-_stress_old[_qp].row(_component)) / _dt * _grad_test[_i][_qp];

}

Real
LagrangianStiffPropDampingImplicit::computeQpJacobian()
{

  if ( _t < _dt)
    return 0.0;
  else
    return TotalLagrangianStressDivergence::computeQpJacobian() * (_q / _dt);
}

Real
LagrangianStiffPropDampingImplicit::computeQpOffDiagJacobian(unsigned int jvar)
{

  if ( _t < _dt)
    return 0.0;
  else
    return TotalLagrangianStressDivergence::computeQpOffDiagJacobian(jvar) * (_q / _dt);
}
