//* This file implements the AD version of enhanced history energy formulation
//* for three-field (u, w, p) poroelastodynamics with phase-field damage
//*
//* H+ = max_{t∈[0,T]} ⟨ 2ψ_o^{e+} + Γp² - ΓM²(∇·w)² + (1-φ_o)ρ^f τ_t/φ² |ẇ|² - (1-φ_o)(ρ^f-ρ^s)|u̇|² ⟩^+

#include "ElkADThreeFieldHistoryEnergyEnhanced.h"

registerMooseObject("farmsApp", ElkADThreeFieldHistoryEnergyEnhanced);

InputParameters
ElkADThreeFieldHistoryEnergyEnhanced::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription(
      "AD version: Computes the enhanced history energy H+ for three-field "
      "poroelastodynamics formulation. Includes pressure, fluid divergence, "
      "and kinetic energy terms for phase-field fracture coupling.");

  // Strain energy
  params.addRequiredParam<MaterialPropertyName>(
      "psie_active", "Name of the positive elastic energy density material property");

  // Pressure variable
  params.addRequiredCoupledVar("pore_pressure", "Pore pressure variable");

  // Solid velocity components
  params.addRequiredCoupledVar("solid_velocity_x", "Solid velocity x-component");
  params.addRequiredCoupledVar("solid_velocity_y", "Solid velocity y-component");
  params.addCoupledVar("solid_velocity_z", "Solid velocity z-component (for 3D)");

  // Fluid velocity components
  params.addRequiredCoupledVar("fluid_velocity_x", "Fluid velocity x-component");
  params.addRequiredCoupledVar("fluid_velocity_y", "Fluid velocity y-component");
  params.addCoupledVar("fluid_velocity_z", "Fluid velocity z-component (for 3D)");

  // Fluid displacement (for divergence calculation)
  params.addRequiredCoupledVar("fluid_disp_x", "Fluid displacement x-component");
  params.addRequiredCoupledVar("fluid_disp_y", "Fluid displacement y-component");
  params.addCoupledVar("fluid_disp_z", "Fluid displacement z-component (for 3D)");

  // Material properties
  params.addRequiredParam<Real>("initial_porosity", "Initial porosity φ_o");
  params.addRequiredParam<Real>("fluid_bulk_modulus", "Fluid bulk modulus K_f");
  params.addRequiredParam<Real>("grain_bulk_modulus", "Grain bulk modulus K_s");
  params.addRequiredParam<MaterialPropertyName>("bulk_modulus", "Solid bulk modulus K");
  params.addParam<MaterialPropertyName>(
      "biot_modulus", "biot_modulus", "Biot modulus M material property name");
  params.addParam<MaterialPropertyName>("porosity", "porosity", "Porosity material property name");

  // Density and tortosity
  params.addRequiredParam<Real>("fluid_density", "Fluid density ρ^f");
  params.addRequiredParam<Real>("solid_density", "Solid density ρ^s");
  params.addRequiredParam<Real>("tortosity", "Tortosity τ_t");

  // Output property name
  params.addParam<MaterialPropertyName>("psie_active_enhanced",
                                        "psie_active_enhanced",
                                        "Name of the enhanced history energy output property");

  // Formulation options
  params.addParam<bool>("use_two_field_formulation",
                        false,
                        "If true, use two-field H+ formulation (only 2*psi_e+ + Gamma*p^2). "
                        "This disables all additional three-field terms.");
  params.addParam<bool>("include_fluid_divergence_term",
                        true,
                        "Include the fluid divergence term: -Gamma*M^2*(div w)^2");
  params.addParam<bool>("include_fluid_kinetic_term",
                        true,
                        "Include the fluid kinetic term: (1-phi_o)*rho_f*tau/phi^2 * |w_dot|^2");
  params.addParam<bool>("include_density_diff_kinetic_term",
                        true,
                        "Include the density difference kinetic term: -(1-phi_o)*(rho_f-rho_s)*|u_dot|^2");

  return params;
}

ElkADThreeFieldHistoryEnergyEnhanced::ElkADThreeFieldHistoryEnergyEnhanced(
    const InputParameters & parameters)
  : ADMaterial(parameters),
    _psie_active(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("psie_active"))),
    _psie_inactive(getADMaterialProperty<Real>("psie_inactive")),
    _g(getADMaterialProperty<Real>("g")),
    _pp(adCoupledValue("pore_pressure")),
    _vel_x(adCoupledValue("solid_velocity_x")),
    _vel_y(adCoupledValue("solid_velocity_y")),
    _vel_z(isCoupled("solid_velocity_z") ? &adCoupledValue("solid_velocity_z") : nullptr),
    _vf_x(adCoupledValue("fluid_velocity_x")),
    _vf_y(adCoupledValue("fluid_velocity_y")),
    _vf_z(isCoupled("fluid_velocity_z") ? &adCoupledValue("fluid_velocity_z") : nullptr),
    _grad_wf_x(adCoupledGradient("fluid_disp_x")),
    _grad_wf_y(adCoupledGradient("fluid_disp_y")),
    _grad_wf_z(isCoupled("fluid_disp_z") ? &adCoupledGradient("fluid_disp_z") : nullptr),
    _initial_porosity(getParam<Real>("initial_porosity")),
    _fluid_bulk_modulus(getParam<Real>("fluid_bulk_modulus")),
    _grain_bulk_modulus(getParam<Real>("grain_bulk_modulus")),
    _K(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("bulk_modulus"))),
    _biot_modulus(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("biot_modulus"))),
    _porosity(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("porosity"))),
    _tortosity(getParam<Real>("tortosity")),
    _fluid_density(getParam<Real>("fluid_density")),
    _solid_density(getParam<Real>("solid_density")),
    _mesh_dimension(_mesh.dimension()),
    _use_two_field_formulation(getParam<bool>("use_two_field_formulation")),
    _include_fluid_divergence_term(getParam<bool>("include_fluid_divergence_term")),
    _include_fluid_kinetic_term(getParam<bool>("include_fluid_kinetic_term")),
    _include_density_diff_kinetic_term(getParam<bool>("include_density_diff_kinetic_term")),
    _psie_active_enhanced(
        declareADProperty<Real>(getParam<MaterialPropertyName>("psie_active_enhanced"))),
    _psie_enhanced(declareADProperty<Real>("psie_enhanced")),
    _psie_active_enhanced_old(
        getMaterialPropertyOld<Real>(getParam<MaterialPropertyName>("psie_active_enhanced"))),
    _psie_enhanced_old(getMaterialPropertyOld<Real>("psie_enhanced"))
{
  // Validate parameters
  if (_initial_porosity <= 0.0 || _initial_porosity >= 1.0)
    paramError("initial_porosity", "Initial porosity must be between 0 and 1");
  if (_fluid_bulk_modulus <= 0.0)
    paramError("fluid_bulk_modulus", "Fluid bulk modulus must be positive");
  if (_grain_bulk_modulus <= 0.0)
    paramError("grain_bulk_modulus", "Grain bulk modulus must be positive");
  if (_tortosity <= 0.0)
    paramError("tortosity", "Tortosity must be positive");
}

void
ElkADThreeFieldHistoryEnergyEnhanced::initQpStatefulProperties()
{
  _psie_active_enhanced[_qp] = 0.0;
  _psie_enhanced[_qp] = 0.0;
}

void
ElkADThreeFieldHistoryEnergyEnhanced::computeQpProperties()
{
  // Compute Gamma coefficient: Γ = (φ_o - 1)(1/K_f - 1/K_s) - K/K_s²
  const ADReal Gamma = computeGamma();

  // Base terms (always included): 2ψ_o^{e+} + Γp²
  ADReal H_current = 2.0 * _psie_active[_qp];        // Elastic energy term
                     //+ Gamma * _pp[_qp] * _pp[_qp]; // Pressure term

  // Additional three-field terms (conditionally included)
  if (!_use_two_field_formulation)
  {
    // Fluid divergence term: -ΓM²(∇·w)²
    if (_include_fluid_divergence_term)
    {
      const ADReal M = _biot_modulus[_qp];
      const ADReal div_w = computeFluidDivergence();
      H_current -= Gamma * M * M * div_w * div_w;
    }

    // Fluid kinetic term: (1-φ_o)ρ^f τ_t/φ² |ẇ|²
    if (_include_fluid_kinetic_term)
    {
      const ADReal phi = _porosity[_qp];
      const ADReal fluid_kinetic_coeff =
          (1.0 - _initial_porosity) * _fluid_density * _tortosity / (phi * phi);
      const ADReal vf_squared = computeFluidVelocitySquared();
      H_current += fluid_kinetic_coeff * vf_squared;
    }

    // Density difference kinetic term: -(1-φ_o)(ρ^f - ρ^s)|u̇|²
    if (_include_density_diff_kinetic_term)
    {
      const Real density_diff_coeff =
          (1.0 - _initial_porosity) * (_fluid_density - _solid_density);
      const ADReal vs_squared = computeSolidVelocitySquared();
      H_current -= density_diff_coeff * vs_squared;
    }
  }

  // Apply Macaulay bracket: ⟨x⟩^+ = max(0, x)
  const ADReal H_current_positive = std::max(ADReal(0.0), H_current);

  // Take maximum over time (irreversibility condition)
  // H+ = max_t ( H_current_positive )
  _psie_active_enhanced[_qp] = std::max(H_current_positive, ADReal(_psie_active_enhanced_old[_qp]));

  // Compute full enhanced strain energy
  _psie_enhanced[_qp] = _g[_qp] * _psie_active_enhanced[_qp] + _psie_inactive[_qp];
}

ADReal
ElkADThreeFieldHistoryEnergyEnhanced::computeGamma() const
{
  // Γ = (φ_o - 1)(1/K_f - 1/K_s) - K/K_s²
  const Real phi_o_minus_1 = _initial_porosity - 1.0;
  const Real inv_Kf = 1.0 / _fluid_bulk_modulus;
  const Real inv_Ks = 1.0 / _grain_bulk_modulus;
  const ADReal K_current = _K[_qp];
  const Real Ks_squared = _grain_bulk_modulus * _grain_bulk_modulus;

  return phi_o_minus_1 * (inv_Kf - inv_Ks) - K_current / Ks_squared;
}

ADReal
ElkADThreeFieldHistoryEnergyEnhanced::computeFluidDivergence() const
{
  // ∇·w = ∂w_x/∂x + ∂w_y/∂y + ∂w_z/∂z
  ADReal div_w = _grad_wf_x[_qp](0) + _grad_wf_y[_qp](1);

  if (_mesh_dimension == 3 && _grad_wf_z)
    div_w += (*_grad_wf_z)[_qp](2);

  return div_w;
}

ADReal
ElkADThreeFieldHistoryEnergyEnhanced::computeFluidVelocitySquared() const
{
  // |ẇ|² = ẇ_x² + ẇ_y² + ẇ_z²
  ADReal vf_squared = _vf_x[_qp] * _vf_x[_qp] + _vf_y[_qp] * _vf_y[_qp];

  if (_mesh_dimension == 3 && _vf_z)
    vf_squared += (*_vf_z)[_qp] * (*_vf_z)[_qp];

  return vf_squared;
}

ADReal
ElkADThreeFieldHistoryEnergyEnhanced::computeSolidVelocitySquared() const
{
  // |u̇|² = u̇_x² + u̇_y² + u̇_z²
  ADReal vs_squared = _vel_x[_qp] * _vel_x[_qp] + _vel_y[_qp] * _vel_y[_qp];

  if (_mesh_dimension == 3 && _vel_z)
    vs_squared += (*_vel_z)[_qp] * (*_vel_z)[_qp];

  return vs_squared;
}
