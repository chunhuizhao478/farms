#pragma once

#include "Kernel.h"

/**
 * This kernel implements the damage evolution forcing term:
 * $(1-B) \psi f(\alpha, xi)$
 * where $f(\alpha, xi)$ is the conditinal forcing term
 * f(\alpha, xi) = Cd * I2 * (xi - xi_0) if xi >= xi_0
 * f(\alpha, xi) = C1 * exp(alpha/C2) * I2 * (xi - xi_0) if xi < xi_0
 */

class BreakageEvolutionConditionalForcing : public Kernel
{
public:
  static InputParameters validParams();

  BreakageEvolutionConditionalForcing(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  
  //implement the conditional forcing term f(\alpha, xi)
  virtual Real computedamageevolutionforcingterm();
  //implement the derivative of the conditional forcing term f(\alpha, xi) with respect to alpha
  virtual Real computedamageevolutionforcingterm_derivative();  

private:
  const MaterialProperty<Real> & _B_breakagevar; //breakage variable
  const MaterialProperty<Real> & _Cd; //damage rate coefficient
  const MaterialProperty<Real> & _I2; //second elastic strain invariant
  const MaterialProperty<Real> & _xi; //strain invariant ratio
  const MaterialProperty<Real> & _xi_0; //strain invariant ratio: onset of damage evolution
  const MaterialProperty<Real> & _C1; //coefficient of healing for damage evolution
  const MaterialProperty<Real> & _C2; //coefficient for exponential term in damage evolution
};