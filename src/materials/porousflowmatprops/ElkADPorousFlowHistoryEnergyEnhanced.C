//* This file implements the AD version of enhanced history energy formulation for phase field damage
//* with poroelastic property evolution as described in the CMAME paper Appendix A

#include "ElkADPorousFlowHistoryEnergyEnhanced.h"

registerMooseObject("farmsApp", ElkADPorousFlowHistoryEnergyEnhanced);

InputParameters
ElkADPorousFlowHistoryEnergyEnhanced::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription(
      "AD version: Computes the enhanced history energy H+ that includes pressure-dependent terms "
      "for poroelastic coupling in phase field fracture. Implements equation (A.12) "
      "from the CMAME paper Appendix A.");

  params.addRequiredParam<MaterialPropertyName>(
      "psie_active", "Name of the positive elastic energy density material property");
  params.addRequiredCoupledVar("pore_pressure", "Pore pressure variable");
  params.addRequiredParam<Real>("initial_porosity", "Initial porosity φ_o");
  params.addRequiredParam<Real>("fluid_bulk_modulus", "Fluid bulk modulus K_f");
  params.addRequiredParam<Real>("grain_bulk_modulus", "Grain bulk modulus K_s");
  params.addRequiredParam<MaterialPropertyName>("bulk_modulus", "Solid bulk modulus K");
  params.addParam<MaterialPropertyName>("psie_active_enhanced",
                                        "psie_active_enhanced",
                                        "Name of the enhanced history energy material property (output)");

  return params;
}

ElkADPorousFlowHistoryEnergyEnhanced::ElkADPorousFlowHistoryEnergyEnhanced(
    const InputParameters & parameters)
  : ADMaterial(parameters),
    _psie_active(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("psie_active"))),
    _psie_inactive(getADMaterialProperty<Real>("psie_inactive")),
    _g(getADMaterialProperty<Real>("g")),
    _pp(adCoupledValue("pore_pressure")),
    _initial_porosity(getParam<Real>("initial_porosity")),
    _fluid_bulk_modulus(getParam<Real>("fluid_bulk_modulus")),
    _grain_bulk_modulus(getParam<Real>("grain_bulk_modulus")),
    _K(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("bulk_modulus"))),
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
}

void
ElkADPorousFlowHistoryEnergyEnhanced::initQpStatefulProperties()
{
  _psie_active_enhanced[_qp] = 0.0;
  _psie_enhanced[_qp] = 0.0;
}

void
ElkADPorousFlowHistoryEnergyEnhanced::computeQpProperties()
{
  // Compute the pressure-dependent coefficient
  // coeff = (φ_o - 1)*(1/K_f - 1/K_s) - K/K_s²
  const ADReal coeff = computePressureCoefficient();

  // Compute the enhanced history energy before applying Macaulay bracket
  // H_current = 2*ψ_e^{e+} + p²*coeff
  const ADReal H_current = 2.0 * _psie_active[_qp] + _pp[_qp] * _pp[_qp] * coeff;

  // Apply Macaulay bracket: <x>^+ = max(0, x)
  const ADReal H_current_positive = std::max(ADReal(0.0), H_current);

  // Take maximum over time (irreversibility condition)
  // H+ = max_t ( H_current_positive )
  _psie_active_enhanced[_qp] = std::max(H_current_positive, ADReal(_psie_active_enhanced_old[_qp]));
  _psie_enhanced[_qp] = _g[_qp] * _psie_active_enhanced[_qp] + _psie_inactive[_qp];
}

ADReal
ElkADPorousFlowHistoryEnergyEnhanced::computePressureCoefficient() const
{
  // Compute: (φ_o - 1)*(1/K_f - 1/K_s) - K/K_s²
  //
  // This coefficient appears in equation (A.12) of the paper's Appendix A.
  // It represents the contribution of poroelastic property evolution to the history energy.

  const Real phi_o_minus_1 = _initial_porosity - 1.0;
  const Real inv_Kf = 1.0 / _fluid_bulk_modulus;
  const Real inv_Ks = 1.0 / _grain_bulk_modulus;
  const ADReal K_current = _K[_qp]; // Current solid bulk modulus at this quadrature point
  const Real Ks_squared = _grain_bulk_modulus * _grain_bulk_modulus;

  const ADReal coeff = phi_o_minus_1 * (inv_Kf - inv_Ks) - K_current / Ks_squared;

  return coeff;
}
