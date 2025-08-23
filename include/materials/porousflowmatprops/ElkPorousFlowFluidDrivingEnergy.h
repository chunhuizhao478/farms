//* This file is part of the FARMS application
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "PorousFlowMaterialVectorBase.h"

/**
 * ElkPorousFlowFluidDrivingEnergy
 *
 * Computes the fluid driving energy density per unit volume for a fully saturated
 * single-phase porous medium, following
 *   psi_f = 0.5 * M * [ (tr(eps))^2 - 2 * alpha * theta * tr(eps) + theta^2 ]
 * where theta = p / M + alpha * tr(eps),
 * M is the Biot modulus, alpha is the Biot coefficient, p is pore pressure,
 * and tr(eps) is the volumetric strain.
 *
 * Exposes AD material property:
 *   - fluid_driving_energy_density
 */
class ElkPorousFlowFluidDrivingEnergy : public PorousFlowMaterialVectorBase
{
public:
  static InputParameters validParams();

  ElkPorousFlowFluidDrivingEnergy(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  // Input parameters/properties
  const Real _biot_coefficient;                   // alpha
  const MaterialProperty<Real> & _M;              // PorousFlow_constant_biot_modulus_qp
  const MaterialProperty<Real> & _eps_v;          // PorousFlow_total_volumetric_strain_qp
  const MaterialProperty<std::vector<Real>> & _p; // PorousFlow_porepressure_qp (per-phase vector)

  // Output non-AD property (use ElementIntegralMaterialProperty to integrate)
  MaterialProperty<Real> & _psi_f;
};
