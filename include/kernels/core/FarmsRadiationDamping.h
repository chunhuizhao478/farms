#pragma once

#include "TimeDerivative.h"

/**
 * Domain-wide radiation damping kernel.
 *
 * Adds a damping force proportional to the time derivative of the
 * displacement variable:  eta * u_dot.  The damping coefficient can
 * be supplied either as a constant parameter or as a material
 * property.  An optional flag property can be provided to enable or
 * disable the damping on a per-quadrature-point basis.
 */
class FarmsRadiationDamping : public TimeDerivative
{
public:
  static InputParameters validParams();

  FarmsRadiationDamping(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Returns the damping coefficient at the current qp
  Real dampingCoefficient() const;

  /// Constant damping coefficient (used when no material property is provided)
  const Real _eta_constant;
  /// Optional spatially varying damping coefficient
  const MaterialProperty<Real> * const _eta_property;
  /// Optional flag (e.g., from ForceDampingFlag) to toggle damping
  const MaterialProperty<Real> * const _flag_property;
};
