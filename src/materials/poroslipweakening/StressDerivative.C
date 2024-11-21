
#include "StressDerivative.h"


registerMooseObject("farmsApp", StressDerivative);

InputParameters
StressDerivative::validParams()
{
  InputParameters params = StressDivergenceTensors::validParams();
  params.addClassDescription("Reaction force derivative.");
  return params;
}

StressDerivative::StressDerivative(const InputParameters & parameters)
  : StressDivergenceTensors(parameters),
    _dreaction(getMaterialPropertyByName<RankTwoTensor>(_base_name + "_dreaction"))
{
}


Real
StressDerivative::computeDR_diag()

{
return StressDivergenceTensors::computeQpJacobian();
}

Real
StressDerivative::computeDR_offdiag(unsigned int jvar)
{
 return StressDivergenceTensors::computeQpOffDiagJacobian(jvar);
}

