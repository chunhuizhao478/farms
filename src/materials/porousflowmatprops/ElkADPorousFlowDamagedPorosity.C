//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkADPorousFlowDamagedPorosity.h"

#include "MooseError.h"

#include <algorithm>
#include <cmath>

registerMooseObject("farmsApp", ElkADPorousFlowDamagedPorosity);

InputParameters
ElkADPorousFlowDamagedPorosity::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addRequiredCoupledVar("phase_field", "Damage (phase-field) variable, d.");
  params.addRequiredRangeCheckedParam<Real>(
      "initial_porosity", "initial_porosity>=0 & initial_porosity<=1", "Undamaged porosity phi_0");
  params.addParam<Real>(
      "porosity_lower_bound",
      0.0,
      "Lower clamp applied to the updated porosity (intact-granite anchor; see "
      "test/tests/materials/damaged_porosity/POROSITY_BOUNDS_REFERENCE.md).");
  params.addParam<Real>(
      "porosity_upper_bound",
      0.999,
      "Upper clamp applied to the updated porosity; the default 0.999 regularizes the "
      "fully-damaged phi->1 void limit (see "
      "test/tests/materials/damaged_porosity/POROSITY_BOUNDS_REFERENCE.md).");
  params.addClassDescription("AD material that computes damage-dependent porosity phi(d) = phi_0 + (1 - phi_0)[1 - (1 - d)^2]");
  return params;
}

ElkADPorousFlowDamagedPorosity::ElkADPorousFlowDamagedPorosity(const InputParameters & parameters)
  : ADMaterial(parameters),
    _damage(adCoupledValue("phase_field")),
    _initial_porosity(getParam<Real>("initial_porosity")),
    _porosity_lower_bound(getParam<Real>("porosity_lower_bound")),
    _porosity_upper_bound(getParam<Real>("porosity_upper_bound")),
    _porosity_damaged(declareADProperty<Real>(isNodal() ? "PorousFlow_porosity_nodal_damaged"
                                                        : "PorousFlow_porosity_qp_damaged"))
{
  if (_porosity_lower_bound > _porosity_upper_bound)
    mooseError("porosity_lower_bound must not exceed porosity_upper_bound in ",
               name(),
               ".");
}

void
ElkADPorousFlowDamagedPorosity::computeQpProperties()
{
  // Clamp damage to [0, 1]
  ADReal dmg = std::min(std::max(_damage[_qp], ADReal(0.0)), ADReal(1.0));

  // Degradation function g(d) = (1 - d)^2 (valid for AT1 or AT2 model)
  ADReal g = pow(1.0 - dmg, 2);

  // Damaged porosity: phi(d) = phi_0 + (1 - phi_0) * (1 - g)
  ADReal phi_raw = _initial_porosity + (1.0 - _initial_porosity) * (1.0 - g);

  // Clamp to bounds
  _porosity_damaged[_qp] = std::min(std::max(phi_raw, ADReal(_porosity_lower_bound)),
                                    ADReal(_porosity_upper_bound));
}
