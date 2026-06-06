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
  params.addParam<MooseEnum>(
      "porosity_update_model",
      MooseEnum("damage strain", "damage"),
      "Porosity update law. 'damage' (default): phi = phi_0 + (1-phi_0)[1-(1-d)^2], the "
      "existing damage-driven bounded-maximum-porosity model. 'strain': phi = phi_0 + eps_1, "
      "the strain-based update of Liu et al. (2024, CMAME 429:117165) eq. (40), where eps_1 "
      "is the maximum principal strain of `strain_property`.");
  params.addParam<MaterialPropertyName>(
      "strain_property",
      "mechanical_strain",
      "Rank-two kinematic strain tensor whose maximum principal value drives the strain-based "
      "porosity update. Only used when porosity_update_model = strain. Defaults to the "
      "'mechanical_strain' property declared by ComputeSmallStrain.");
  params.addClassDescription(
      "Computes the damaged porosity property PorousFlow_porosity_{qp,nodal}_damaged using "
      "either the damage-driven bounded-maximum law phi = phi0 + (1-phi0)[1-(1-d)^2] or the "
      "strain-based law phi = phi0 + eps_1 (Liu et al. 2024 CMAME eq. 40), clamped to "
      "[porosity_lower_bound, porosity_upper_bound].");
  return params;
}

ElkPorousFlowDamagedPorosity::ElkPorousFlowDamagedPorosity(const InputParameters & parameters)
  : Material(parameters),
    _porosity_update_model(getParam<MooseEnum>("porosity_update_model") == "strain"
                               ? PorosityUpdateModel::STRAIN
                               : PorosityUpdateModel::DAMAGE),
    _damage(coupledValue("phase_field")),
    _initial_porosity(getParam<Real>("initial_porosity")),
    _porosity_lower_bound(getParam<Real>("porosity_lower_bound")),
    _porosity_upper_bound(getParam<Real>("porosity_upper_bound")),
    _mechanical_strain(nullptr),
    _porosity_damaged(declareProperty<Real>(isNodal() ? "PorousFlow_porosity_nodal_damaged"
                                                      : "PorousFlow_porosity_qp_damaged"))
{
  if (_porosity_lower_bound > _porosity_upper_bound)
    mooseError("porosity_lower_bound must not exceed porosity_upper_bound in ",
               name(),
               ".");

  if (_porosity_update_model == PorosityUpdateModel::STRAIN)
  {
    // mechanical_strain is a quadrature-point property; it cannot be evaluated for
    // the nodal porosity instance, so the strain-based law is qp-only.
    if (isNodal())
      mooseError("porosity_update_model = strain is only available at quadrature points; the "
                 "qp strain tensor is unavailable for the nodal porosity property in ",
                 name(),
                 ". Use porosity_update_model = damage for nodal porosity.");

    // Bind the strain only in this branch so the default (damage) model does not create a
    // dependency on a strain material that need not be present.
    _mechanical_strain = &getMaterialProperty<RankTwoTensor>("strain_property");
  }
}

void
ElkPorousFlowDamagedPorosity::computeQpProperties()
{
  Real phi_raw;
  if (_porosity_update_model == PorosityUpdateModel::STRAIN)
  {
    // Liu et al. (2024) CMAME 429:117165 eq. (40): phi = phi_0 + eps_1, with eps_1 the maximum
    // (most-tensile) principal strain. symmetricEigenvalues returns the eigenvalues in ascending
    // order (LAPACK dsyev), so the last entry is the maximum principal strain.
    std::vector<Real> eigvals(3); // RankTwoTensor is 3x3; symmetricEigenvalues returns 3, ascending
    (*_mechanical_strain)[_qp].symmetricEigenvalues(eigvals);
    const Real eps1 = eigvals.back(); // ascending order -> last is the maximum principal strain
    phi_raw = _initial_porosity + eps1;
  }
  else
  {
    const Real dmg = std::clamp(_damage[_qp], 0.0, 1.0);
    const Real g = std::pow(1.0 - dmg, 2); // only valid for AT1 or AT2 model
    phi_raw = _initial_porosity + (1.0 - _initial_porosity) * (1.0 - g);
  }
  _porosity_damaged[_qp] = std::clamp(phi_raw, _porosity_lower_bound, _porosity_upper_bound);
}
