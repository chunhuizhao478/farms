#include "StabilizationMaterial.h"

registerMooseObject("farmsApp", StabilizationMaterial);

InputParameters
StabilizationMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Computes PSPG stabilization parameter tau = C * h^2 / (2*G)");
  
  params.addParam<Real>("stabilization_coeff", 0.5, 
                        "Stabilization coefficient C (typical range: 0.1 to 2.0)");
  
  // Add parameter to specify which shear modulus property to use
  params.addParam<MaterialPropertyName>("shear_modulus", "shear_modulus",
                                         "Name of shear modulus material property");
  
  return params;
}

StabilizationMaterial::StabilizationMaterial(const InputParameters & parameters)
  : Material(parameters),
    _tau_pspg(declareProperty<Real>("tau_pspg")),
    _shear_modulus(getMaterialProperty<Real>(getParam<MaterialPropertyName>("shear_modulus"))),
    _coeff(getParam<Real>("stabilization_coeff"))
{
}

void
StabilizationMaterial::computeQpProperties()
{
  // Compute characteristic element length h
  Real h = std::pow(_current_elem->volume(), 1.0 / _mesh.dimension());
  
  // Safety check: ensure shear modulus is positive
  Real G = std::max(_shear_modulus[_qp], 1e-10);
  
  // PSPG stabilization parameter: tau = C * h^2 / (2*G)
  _tau_pspg[_qp] = _coeff * h * h / (2.0 * G);
}