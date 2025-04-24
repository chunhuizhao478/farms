#pragma once

#include "Kernel.h"

/**
 * This kernel implements the damage evolution forcing term:
 * $(1-B) \psi f(\alpha, xi)$
 * where $f(\alpha, xi)$ is the conditinal forcing term
 * f(B, xi) = Cd * I2 * (xi - xi_0) if xi >= xi_0
 * f(B, xi) = C1 * exp(alpha/C2) * I2 * (xi - xi_0) if xi < xi_0
 */

class BreakageEvolutionConditionalForcing : public Kernel
{
public:
  static InputParameters validParams();

  BreakageEvolutionConditionalForcing(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  
  //implement the conditional forcing term f(B, xi)
  virtual Real computebreakageevolutionforcingterm();
  //implement the derivative of the conditional forcing term f(B, xi) with respect to alpha
  virtual Real computebreakageevolutionforcingterm_derivative();  
  //implement the derivative of the conditional forcing term f(B, xi) with respect to alpha
  virtual Real computebreakageevolutionforcingterm_derivative_walpha();
  // Function for alpha_func_root1
  virtual Real alphacr_root1(Real xi);
  // Function for alpha_func_root2
  virtual Real alphacr_root2(Real xi);

private:
  unsigned int _alpha_var;
  const VariableValue & _alpha; //damage variable
  const MaterialProperty<Real> & _Cd; //damage rate coefficient
  const MaterialProperty<Real> & _CdCb_multiplier; //multiplier between Cd and Cb
  const MaterialProperty<Real> & _I2; //second elastic strain invariant
  const MaterialProperty<Real> & _xi; //strain invariant ratio
  const MaterialProperty<Real> & _xi_0; //strain invariant ratio: onset of damage evolution
  const MaterialProperty<Real> & _xi_1; //strain invariant ratio: onset of damage evolution
  const MaterialProperty<Real> & _xi_d; //strain invariant ratio: onset of damage evolution
  const MaterialProperty<Real> & _xi_min; //minimum strain invariant ratio
  const MaterialProperty<Real> & _xi_max; //maximum strain invariant ratio
  const MaterialProperty<Real> & _C1; //coefficient of healing for damage evolution
  const MaterialProperty<Real> & _C2; //coefficient for exponential term in damage evolution
  const MaterialProperty<Real> & _lambda_o; //initial lambda value
  const MaterialProperty<Real> & _shear_modulus_o; //initial shear modulus value
  const MaterialProperty<Real> & _CBH_constant; //constant CBH value
  const MaterialProperty<Real> & _beta_width; //coefficient gives width of transitional region
  const MaterialProperty<Real> & _gamma_damaged_r; //coefficient of healing for breakage evolution
};