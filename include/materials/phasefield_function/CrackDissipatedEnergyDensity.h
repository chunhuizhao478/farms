#pragma once

#include "ADMaterial.h"

/**
 * CrackDissipatedEnergyDensity
 *
 * Computes the dissipated energy density for a phase-field fracture model in the
 * pure solid case (no fluid), following the common crack surface density form:
 *   gamma(d, grad d) = (1/c0) * ( alpha(d)/l + l * |grad d|^2 )
 * The total dissipated energy is then:  Psi^d = ∫ Gc * gamma dV
 *
 * This is equivalent to the alternative representation using alpha'(d) and the
 * Laplacian of d via integration by parts:
 *   gamma = (1/c0) * (1/l * alpha'(d) - 2 l * Δd),
 * but this class uses the grad-based expression to avoid requiring second
 * derivatives of d in the finite element discretization.
 *
 * Inputs (as MaterialProperties or coupled variables):
 *  - phase_field (variable): damage variable d
 *  - alpha (MaterialProperty<Real>): geometric crack function α(d)
 *  - Gc (MaterialProperty<Real>): fracture toughness
 *  - l (MaterialProperty<Real>): regularization length scale
 *  - c0 (MaterialProperty<Real>): normalization constant corresponding to α
 *
 * Output:
 *  - property_name (ADMaterialProperty<Real>, default: "dissipated_energy_density")
 *    equal to Gc/c0 * ( alpha/l + l * |grad d|^2 ). Units: energy per volume
 */
class CrackDissipatedEnergyDensity : public ADMaterial
{
public:
  static InputParameters validParams();
  CrackDissipatedEnergyDensity(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  // Coupled damage variable and its gradient
  const ADVariableValue & _d;
  const ADVariableGradient & _grad_d;

  // Required material properties
  const ADMaterialProperty<Real> & _alpha; // α(d)
  const ADMaterialProperty<Real> & _Gc;    // fracture toughness
  const ADMaterialProperty<Real> & _l;     // length scale
  const ADMaterialProperty<Real> & _c0;    // normalization constant

  // Output property: dissipated energy density per unit volume
  ADMaterialProperty<Real> & _psi_diss;
};
