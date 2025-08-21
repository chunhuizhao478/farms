#include "CrackDissipatedEnergyDensity.h"
#include "MooseVariable.h"
#include "MooseTypes.h"
#include "libmesh/utility.h"

registerMooseObject("farmsApp", CrackDissipatedEnergyDensity);

InputParameters CrackDissipatedEnergyDensity::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription(
      "Computes dissipated energy density Psi^d = Gc/c0 * ( alpha/l + l*|grad d|^2 ) for a pure "
      "solid phase-field fracture model.");

  params.addRequiredCoupledVar("phase_field", "Damage (phase-field) variable d.");
  params.addRequiredParam<MaterialPropertyName>("alpha_name", "Geometric crack function α(d)");
  params.addRequiredParam<MaterialPropertyName>("Gc_name","Fracture toughness Gc");
  params.addRequiredParam<MaterialPropertyName>("l_name","Regularization length l");
  params.addRequiredParam<MaterialPropertyName>("c0_name","Normalization constant c0");
  params.addParam<MaterialPropertyName>(
      "property_name", "dissipated_energy_density", "Output property name for Psi^d");
  return params;
}

CrackDissipatedEnergyDensity::CrackDissipatedEnergyDensity(const InputParameters & parameters)
  : ADMaterial(parameters),
    _d(adCoupledValue("phase_field")),
    _grad_d(adCoupledGradient("phase_field")),
    _alpha(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("alpha_name"))),
    _Gc(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("Gc_name"))),
    _l(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("l_name"))),
    _c0(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("c0_name"))),
    _psi_diss(declareADProperty<Real>(getParam<MaterialPropertyName>("property_name")))
{
}

void CrackDissipatedEnergyDensity::computeQpProperties()
{
  // gamma(d, grad d) = (1/c0) * ( alpha/l + l * |grad d|^2 )
  const ADReal l = _l[_qp];
  const ADReal inv_c0 = 1.0 / _c0[_qp];
  const ADReal alpha_over_l = _alpha[_qp] / l;
  const ADReal grad_term = l * (_grad_d[_qp] * _grad_d[_qp]); // dot product
  const ADReal gamma = inv_c0 * (alpha_over_l + grad_term);

  // Psi^d density = Gc * gamma
  _psi_diss[_qp] = _Gc[_qp] * gamma;
}
