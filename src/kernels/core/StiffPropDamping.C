/*
Kernel of Stiffness Proportional Damping 
*/

#include "StiffPropDamping.h"
#include "ElasticityTensorTools.h"

registerMooseObject("farmsApp", StiffPropDamping);

InputParameters
StiffPropDamping::validParams()
{
  InputParameters params = StressDivergenceTensors::validParams();
  params.addClassDescription(
      "Compute Stiffness Proportional Damping Residual");

  params.addRequiredParam<Real>("q","Ratio Factor to assign magnitude of stiffness proportional damping term");

  return params;
}

StiffPropDamping::StiffPropDamping(const InputParameters & parameters)
  : StressDivergenceTensors(parameters),
  
  //Get stress tensor from previous time step (older == previous time step in Explicit Time Integration Scheme)
  _stress_older(getMaterialPropertyOlderByName<RankTwoTensor>(_base_name + "stress")),
  
  //Get stress tensor from current time step
  _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
  
  //Ratio factor 
  _q(getParam<Real>("q"))

{
}

Real
StiffPropDamping::computeQpResidual()
{

  Real residual = _q*(_stress[_qp].row(_component)-_stress_older[_qp].row(_component)) * _grad_test[_i][_qp];

  return residual;
}

Real
StiffPropDamping::computeQpJacobian()
{

  return StressDivergenceTensors::computeQpJacobian()*_q;

}

Real
StiffPropDamping::computeQpOffDiagJacobian(unsigned int jvar)
{

  return 0.0;
  
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    if (jvar == _disp_var[i])
    { 
      return StressDivergenceTensors::computeQpOffDiagJacobian(jvar)*_q;   
    }     
  }
}



