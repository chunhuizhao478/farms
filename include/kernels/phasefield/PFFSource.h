//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "KernelValue.h"
#include "DerivativeMaterialPropertyNameInterface.h"

class PFFSource : public KernelValue, public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  PFFSource(const InputParameters & parameters);

protected:
  virtual Real precomputeQpResidual() override;
  virtual Real precomputeQpJacobian() override;

  /// The derivative of Helmholtz free energy w.r.t. the phase field
  const MaterialProperty<Real> & _dpsi_dd;

  /// The second derivative of Helmholtz free energy w.r.t. the phase field
  const MaterialProperty<Real> & _d2psi_dd2;
};
