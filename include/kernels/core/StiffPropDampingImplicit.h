#pragma once

#include "TotalLagrangianStressDivergence.h"

class StiffPropDampingImplicit : public TotalLagrangianStressDivergence
{
public:
  static InputParameters validParams();

  StiffPropDampingImplicit(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<RankTwoTensor> & _stress_old;
  const MaterialProperty<RankTwoTensor> & _stress;
  Real _q;
  /// An integer corresponding to the direction this kernel acts in
  const unsigned int _component;
};