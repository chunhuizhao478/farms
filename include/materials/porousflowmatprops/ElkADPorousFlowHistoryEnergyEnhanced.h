//* This file implements the AD version of enhanced history energy formulation for phase field damage
//* with poroelastic property evolution as described in the CMAME paper Appendix A

#pragma once

#include "ADMaterial.h"

/**
 * AD version of Material that computes the enhanced history energy H+ including pressure-dependent
 * terms for poroelastic coupling in phase field fracture.
 *
 * From the paper's Appendix A, equation (A.12):
 * H+ = max_t ( (2*ψ_e^{e+} + p²*[(φ_o - 1)*(1/K_f - 1/K_s) - K/K_s²])^+ )
 *
 * Where:
 * - ψ_e^{e+} is the positive elastic energy density
 * - p is pore pressure
 * - φ_o is initial porosity
 * - K_f is fluid bulk modulus
 * - K_s is grain bulk modulus
 * - K is solid bulk modulus
 */
class ElkADPorousFlowHistoryEnergyEnhanced : public ADMaterial
{
public:
  static InputParameters validParams();

  ElkADPorousFlowHistoryEnergyEnhanced(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Positive elastic energy density (from ADSmallDeformationIsotropicElasticityPF)
  const ADMaterialProperty<Real> & _psie_active;
  const ADMaterialProperty<Real> & _psie_inactive;
  const ADMaterialProperty<Real> & _g;

  /// Pore pressure
  const ADVariableValue & _pp;

  /// Initial porosity
  const Real _initial_porosity;

  /// Fluid bulk modulus
  const Real _fluid_bulk_modulus;

  /// Grain bulk modulus
  const Real _grain_bulk_modulus;

  /// Solid bulk modulus (from material property)
  const ADMaterialProperty<Real> & _K;

  /// Enhanced history energy (output)
  ADMaterialProperty<Real> & _psie_active_enhanced;
  ADMaterialProperty<Real> & _psie_enhanced;

  /// Old value of enhanced history energy (for max operation)
  const MaterialProperty<Real> & _psie_active_enhanced_old;
  const MaterialProperty<Real> & _psie_enhanced_old;

  /// Compute pressure coefficient for the additional term
  ADReal computePressureCoefficient() const;
};
