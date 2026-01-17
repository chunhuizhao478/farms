//* This file is part of the FARMS application
//* Rate-dependent phase field fracture kernel
//* Based on Hofacker & Miehe (2012) - IJNME 93:276-301

#pragma once

#include "ADKernelValue.h"

/**
 * ADPFFViscousResistance implements the viscous crack resistance term
 * for rate-dependent phase field fracture.
 *
 * This kernel adds the term: eta * d_dot
 * where:
 *   - eta is the viscosity parameter [Pa*s] or [N*s/m^2]
 *   - d_dot is the time derivative of the phase field
 *
 * Using backward Euler time discretization:
 *   d_dot = (d - d_old) / dt
 *
 * Weak form contribution:
 *   (w, eta * d_dot) = (w, eta/dt * (d - d_old))
 *
 * Reference:
 *   Hofacker, M., & Miehe, C. (2013). A phase field model of dynamic fracture:
 *   Robust field updates for the analysis of complex crack patterns.
 *   Int. J. Numer. Meth. Engng, 93:276-301. Equation (46).
 */
class ADPFFViscousResistance : public ADKernelValue
{
public:
  static InputParameters validParams();

  ADPFFViscousResistance(const InputParameters & parameters);

protected:
  virtual ADReal precomputeQpResidual() override;

  /// The viscosity parameter eta [Pa*s] or [N*s/m^2]
  const ADMaterialProperty<Real> & _eta;

  /// Old value of phase field from previous time step
  const VariableValue & _d_old;
};
