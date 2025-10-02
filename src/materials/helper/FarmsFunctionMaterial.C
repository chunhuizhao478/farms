#include "FarmsFunctionMaterial.h"

#include "Function.h"

registerMooseObject("farmsApp", FarmsFunctionMaterial);

InputParameters
FarmsFunctionMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Copies the value of a Function into a scalar material property.");
  params.addRequiredParam<MaterialPropertyName>("property_name",
                                                "Name of the material property to populate");
  params.addRequiredParam<FunctionName>("function",
                                        "Function that will be evaluated at the quadrature points");
  return params;
}

FarmsFunctionMaterial::FarmsFunctionMaterial(const InputParameters & parameters)
  : Material(parameters),
    _function(getFunction("function")),
    _property(declareProperty<Real>(getParam<MaterialPropertyName>("property_name")))
{
}

void
FarmsFunctionMaterial::initQpStatefulProperties()
{
  // Initialize the properties
  _property[_qp] = _function.value(_t, _q_point[_qp]);
}

void
FarmsFunctionMaterial::computeQpProperties()
{
  _property[_qp] = _function.value(_t, _q_point[_qp]);
}

