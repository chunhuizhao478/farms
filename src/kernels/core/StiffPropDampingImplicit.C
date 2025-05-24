/*
Kernel of Stiffness Proportional Damping 
*/

#include "StiffPropDampingImplicit.h"
#include "ElasticityTensorTools.h"

registerMooseObject("farmsApp", StiffPropDampingImplicit);

InputParameters
StiffPropDampingImplicit::validParams()
{
  InputParameters params = StressDivergenceTensors::validParams();
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

StiffPropDampingImplicit::StiffPropDampingImplicit(const InputParameters & parameters)
  : StressDivergenceTensors(parameters),
  
  //Get stress tensor from previous time step
  _stress_old(getMaterialPropertyOldByName<RankTwoTensor>("stress")),
  
  //Get stress tensor from current time step
  _stress(getMaterialPropertyByName<RankTwoTensor>("stress")),
  
  //Ratio factor 
  _q(getParam<Real>("zeta")),

  //component
  _component(getParam<unsigned int>("component"))

{
}

Real
StiffPropDampingImplicit::computeQpResidual()
{
  if ( _t < _dt)
    return 0.0;
  else
    return _q*(_stress[_qp].row(_component)-_stress_old[_qp].row(_component)) * _grad_test[_i][_qp];

}

Real
StiffPropDampingImplicit::computeQpJacobian()
{

  if ( _t < _dt)
    return 0.0;
  else
    return StressDivergenceTensors::computeQpJacobian() * (_q);
}

Real
StiffPropDampingImplicit::computeQpOffDiagJacobian(unsigned int jvar)
{

  if ( _t < _dt)
    return 0.0;
  else
    return StressDivergenceTensors::computeQpOffDiagJacobian(jvar) * (_q);
}
