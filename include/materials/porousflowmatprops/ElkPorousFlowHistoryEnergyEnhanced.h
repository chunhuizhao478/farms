//* This file implements the enhanced history energy formulation for phase field damage
//* with poroelastic property evolution as described in the CMAME paper Appendix A

#pragma once

#include "Material.h"

/**
 * This Material computes the enhanced history energy H+ that includes pressure-dependent terms
 * for poroelastic coupling in phase field fracture.
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
class ElkPorousFlowHistoryEnergyEnhanced : public Material
{
public:
  static InputParameters validParams();

  ElkPorousFlowHistoryEnergyEnhanced(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Positive elastic energy density (from NDSmallDeformationIsotropicElasticity)
  const MaterialProperty<Real> & _psie_active;
  const MaterialProperty<Real> & _psie_inactive;
  const MaterialProperty<Real> & _g;

  /// Pore pressure
  const VariableValue & _pp;

  /// Initial porosity
  const Real _initial_porosity;

  /// Fluid bulk modulus
  const Real _fluid_bulk_modulus;

  /// Grain bulk modulus
  const Real _grain_bulk_modulus;

  /// Solid bulk modulus (from material property)
  const MaterialProperty<Real> & _K;

  /// Enhanced history energy (output)
  MaterialProperty<Real> & _psie_active_enhanced;
  MaterialProperty<Real> & _psie_enhanced;

  /// Old value of enhanced history energy (for max operation)
  const MaterialProperty<Real> & _psie_active_enhanced_old;
  const MaterialProperty<Real> & _psie_enhanced_old;

  /// Pressure coefficient for the additional term
  Real computePressureCoefficient() const;
};
