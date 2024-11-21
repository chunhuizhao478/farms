
#pragma once

#include "DerivativeFunctionMaterialBase.h"
#include "ADRankTwoTensorForward.h"
#include "ADRankFourTensorForward.h"

/**
 * Material class to compute derivatives of stress
 */
class StressDerivative2 : public DerivativeFunctionMaterialBase
{
public:
  static InputParameters validParams();

  StressDerivative2(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual Real computeDF(unsigned int i_var) override;

  const std::string _base_name;

  /// Stress tensor
  const MaterialProperty<RankTwoTensor> & _stress;
  MaterialProperty<RankTwoTensor> & _dstress;

  ///@{ Elasticity tensor derivatives
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  std::vector<const MaterialProperty<RankFourTensor> *> _delasticity_tensor;

  const MaterialProperty<RankFourTensor> & _Jacobian_mult;

  ///@{ Strain and derivatives
  const MaterialProperty<RankTwoTensor> & _strain;
  std::vector<const MaterialProperty<RankTwoTensor> *> _dstrain;

};