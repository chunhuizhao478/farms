//* This file is part of the FARMS application
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "ElkPorousFlowFluidDrivingEnergy.h"

// no extra includes needed

registerMooseObject("farmsApp", ElkPorousFlowFluidDrivingEnergy);

InputParameters
ElkPorousFlowFluidDrivingEnergy::validParams()
{
  InputParameters params = PorousFlowMaterialVectorBase::validParams();
  params.addClassDescription("Fluid driving energy density psi_f for porous flow HM coupling");
  params.addRequiredRangeCheckedParam<Real>(
      "biot_coefficient",
      "biot_coefficient>=0 & biot_coefficient<=1",
      "Biot coefficient (alpha)");
  // Expect the Biot modulus from PorousFlowConstantBiotModulus or equivalent
  params.set<std::string>("pf_material_type") = "fluid_driving_energy";
  return params;
}

ElkPorousFlowFluidDrivingEnergy::ElkPorousFlowFluidDrivingEnergy(
    const InputParameters & parameters)
  : PorousFlowMaterialVectorBase(parameters),
    _biot_coefficient(getParam<Real>("biot_coefficient")),
    _M(getMaterialProperty<Real>(
        _nodal_material ? "PorousFlow_constant_biot_modulus_nodal"
                        : "PorousFlow_constant_biot_modulus_qp")),
    _eps_v(getMaterialProperty<Real>(
        _nodal_material ? "PorousFlow_total_volumetric_strain_nodal"
                        : "PorousFlow_total_volumetric_strain_qp")),
    _p(getMaterialProperty<std::vector<Real>>(
        _nodal_material ? "PorousFlow_porepressure_nodal" : "PorousFlow_porepressure_qp")),
  _psi_f(declareProperty<Real>("fluid_driving_energy_density"))
{
}

void
ElkPorousFlowFluidDrivingEnergy::computeQpProperties()
{
  // Single-phase assumption: take phase 0 pressure
  const Real p = _p[_qp].empty() ? 0.0 : _p[_qp][0];
  const Real M = _M[_qp];
  const Real tr_eps = _eps_v[_qp];
  const Real alpha = _biot_coefficient;

  // theta = p/M + alpha * tr(eps)
  const Real theta = (M > 0.0 ? p / M : 0.0) + alpha * tr_eps;

  if (p == 0.0)
  {
    std::cout << "Warning: pore pressure is zero at qp " << _qp << std::endl;
  }

  // psi_f = 0.5 * M * [ alpha*alpha*(tr eps)^2 - 2*alpha*theta*tr(eps) + theta^2 ]
  const Real term1 = alpha * alpha * tr_eps * tr_eps;
  const Real term2 = -2.0 * alpha * theta * tr_eps;
  const Real term3 = theta * theta;

  _psi_f[_qp] = 0.5 * M * (term1 + term2 + term3);
}
