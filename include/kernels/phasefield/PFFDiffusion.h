//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Kernel.h"
#include "BaseNameInterface.h"

class PFFDiffusion : public Kernel, public BaseNameInterface
{
public:
  static InputParameters validParams();

  PFFDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// The fracture toughness
  const MaterialProperty<Real> & _Gc;

  /// The normalization constant
  const MaterialProperty<Real> & _c0;

  /// The regularization length
  const MaterialProperty<Real> & _l;
};
