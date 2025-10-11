#pragma once

#include "Kernel.h"

/**
 * Kernel that applies the Sochacki sponge damping contribution to the
 * momentum balance. It introduces a term proportional to the time derivative
 * of the primary variable and the spatially varying attenuation coefficient
 * computed by SochackiSpongeMaterial.
 */
class SochackiSpongeDamping : public Kernel
{
public:
  static InputParameters validParams();
  SochackiSpongeDamping(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

private:
  const VariableValue & _u_dot;
  const MaterialProperty<Real> & _density;
  const MaterialProperty<Real> & _sponge_coeff;
  const Real _damping_scale;
};
