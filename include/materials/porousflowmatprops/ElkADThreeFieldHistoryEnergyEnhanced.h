//* This file implements the AD version of enhanced history energy formulation
//* for three-field (u, w, p) poroelastodynamics with phase-field damage
//*
//* Two-field H+: H+ = max_{t∈[0,T]} ⟨ 2ψ_o^{e+} + Γp² ⟩^+
//*
//* Three-field H+: H+ = max_{t∈[0,T]} ⟨ 2ψ_o^{e+} + Γp² - ΓM²(∇·w)² + (1-φ_o)ρ^f τ_t/φ² |ẇ|² - (1-φ_o)(ρ^f-ρ^s)|u̇|² ⟩^+
//*
//* where: Γ = (φ_o - 1)(1/K_f - 1/K_s) - K/K_s²
//*
//* Use `use_two_field_formulation = true` to disable the additional three-field terms

#pragma once

#include "ADMaterial.h"

/**
 * AD Material that computes the enhanced history energy H+ for three-field
 * poroelastodynamics formulation with phase-field damage coupling.
 *
 * The three-field formulation includes additional terms compared to two-field:
 * - Divergence of fluid displacement: -ΓM²(∇·w)²
 * - Fluid kinetic energy: (1-φ_o)ρ^f τ_t/φ² |ẇ|²
 * - Density difference kinetic: -(1-φ_o)(ρ^f-ρ^s)|u̇|²
 *
 * Set use_two_field_formulation=true to use only the first two terms (2ψ_e+ + Γp²)
 */
class ElkADThreeFieldHistoryEnergyEnhanced : public ADMaterial
{
public:
  static InputParameters validParams();

  ElkADThreeFieldHistoryEnergyEnhanced(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Positive elastic energy density (from elasticity material)
  const ADMaterialProperty<Real> & _psie_active;
  const ADMaterialProperty<Real> & _psie_inactive;
  const ADMaterialProperty<Real> & _g;

  /// Pore pressure
  const ADVariableValue & _pp;

  /// Solid displacement velocities
  const ADVariableValue & _vel_x;
  const ADVariableValue & _vel_y;
  const ADVariableValue * _vel_z;

  /// Fluid displacement velocities
  const ADVariableValue & _vf_x;
  const ADVariableValue & _vf_y;
  const ADVariableValue * _vf_z;

  /// Fluid displacement gradients (for divergence)
  const ADVariableGradient & _grad_wf_x;
  const ADVariableGradient & _grad_wf_y;
  const ADVariableGradient * _grad_wf_z;

  /// Initial porosity
  const Real _initial_porosity;

  /// Fluid bulk modulus
  const Real _fluid_bulk_modulus;

  /// Grain bulk modulus
  const Real _grain_bulk_modulus;

  /// Solid bulk modulus (from material property)
  const ADMaterialProperty<Real> & _K;

  /// Biot modulus M (from material property)
  const ADMaterialProperty<Real> & _biot_modulus;

  /// Porosity (can be damaged)
  const ADMaterialProperty<Real> & _porosity;

  /// Tortosity
  const Real _tortosity;

  /// Fluid density
  const Real _fluid_density;

  /// Solid density
  const Real _solid_density;

  /// Mesh dimension
  const unsigned int _mesh_dimension;

  /// Use two-field formulation (only elastic + pressure terms)
  const bool _use_two_field_formulation;

  /// Include individual three-field terms (for fine-grained control)
  const bool _include_fluid_divergence_term;
  const bool _include_fluid_kinetic_term;
  const bool _include_density_diff_kinetic_term;

  /// Enhanced history energy (output)
  ADMaterialProperty<Real> & _psie_active_enhanced;
  ADMaterialProperty<Real> & _psie_enhanced;

  /// Old value of enhanced history energy (for max operation)
  const MaterialProperty<Real> & _psie_active_enhanced_old;
  const MaterialProperty<Real> & _psie_enhanced_old;

  /// Compute Gamma coefficient: Γ = (φ_o - 1)(1/K_f - 1/K_s) - K/K_s²
  ADReal computeGamma() const;

  /// Compute divergence of fluid displacement: ∇·w
  ADReal computeFluidDivergence() const;

  /// Compute squared magnitude of fluid velocity: |ẇ|²
  ADReal computeFluidVelocitySquared() const;

  /// Compute squared magnitude of solid velocity: |u̇|²
  ADReal computeSolidVelocitySquared() const;
};
