#pragma once

#include "StressDivergenceTensors.h"

class StiffPropDamping : public StressDivergenceTensors
{
public:
  static InputParameters validParams();

  StiffPropDamping(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<RankTwoTensor> & _stress_older;
  const MaterialProperty<RankTwoTensor> & _stress;
  Real _q;

};