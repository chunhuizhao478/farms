#include "ElkPorousFlowDamagedPorosity.h"

#include "MooseError.h"

#include <algorithm>
#include <cmath>

registerMooseObject("farmsApp", ElkPorousFlowDamagedPorosity);

InputParameters
ElkPorousFlowDamagedPorosity::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("phase_field", "Damage (phase-field) variable, c.");
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
  params.addClassDescription("Computes damage-dependent porosity phi(c) = phi0 + (1 - phi0)[1 - (1 - c)^2]");
  return params;
}

ElkPorousFlowDamagedPorosity::ElkPorousFlowDamagedPorosity(const InputParameters & parameters)
  : Material(parameters),
    _damage(coupledValue("phase_field")),
    _initial_porosity(getParam<Real>("initial_porosity")),
    _porosity_lower_bound(getParam<Real>("porosity_lower_bound")),
    _porosity_upper_bound(getParam<Real>("porosity_upper_bound")),
    _porosity_damaged(declareProperty<Real>(isNodal() ? "PorousFlow_porosity_nodal_damaged"
                                                      : "PorousFlow_porosity_qp_damaged"))
{
  if (_porosity_lower_bound > _porosity_upper_bound)
    mooseError("porosity_lower_bound must not exceed porosity_upper_bound in ",
               name(),
               ".");
}

void
ElkPorousFlowDamagedPorosity::computeQpProperties()
{
  const Real dmg = std::clamp(_damage[_qp], 0.0, 1.0);
  const Real g = std::pow(1.0 - dmg, 2); //only valid for AT1 or AT2 model
  const Real phi_raw = _initial_porosity + (1.0 - _initial_porosity) * (1.0 - g);
  _porosity_damaged[_qp] = std::clamp(phi_raw, _porosity_lower_bound, _porosity_upper_bound);
}
