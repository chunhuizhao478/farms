#pragma once
#include "TimeDerivative.h"

/** RhoTimeDerivative
 * rho * dt(u) weak form contribution: integral rho * test * u_dot dOmega.
 * Density material property name provided by parameter rho_name.
 */
class RhoTimeDerivative : public TimeDerivative
{
public:
  static InputParameters validParams();
  RhoTimeDerivative(const InputParameters & p);
protected:
  Real computeQpResidual() override;
  Real computeQpJacobian() override; // simple scaling of TimeDerivative Jacobian
  const MaterialProperty<Real> & _rho;
};
