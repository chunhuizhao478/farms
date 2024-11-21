#pragma once

#include "StressDivergenceTensors.h"

class StressDerivative : public StressDivergenceTensors
{
public:
  static InputParameters validParams();

  StressDerivative(const InputParameters & parameters);

protected:

  virtual Real computeDR_diag();
  virtual Real computeDR_offdiag(unsigned int jvar);

  const MaterialProperty<RankTwoTensor> & _dreaction;


};
