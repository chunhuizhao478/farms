#pragma once

#include "Kernel.h"

/// \brief  Injects the `damage_perturbation` material property
///         as a volumetric source into the α‐equation.
class PerturbationSource : public Kernel
{
public:
  static InputParameters validParams();

  PerturbationSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// This reads the QP‐wise property your Material set
  const MaterialProperty<Real> & _damage_src;
};