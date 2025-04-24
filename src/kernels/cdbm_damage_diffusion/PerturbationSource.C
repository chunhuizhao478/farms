#include "PerturbationSource.h"
#include "MooseVariable.h"

registerMooseObject("farmsApp", PerturbationSource);

InputParameters
PerturbationSource::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<MaterialPropertyName>(
    "damage_source",
    "Name of the material property that holds the damage‐rate source (e.g. damage_perturbation)"
  );
  return params;
}

PerturbationSource::PerturbationSource(const InputParameters & params) :
  Kernel(params),
  _damage_src(getMaterialProperty<Real>(params.get<MaterialPropertyName>("damage_source")))
{}

Real
PerturbationSource::computeQpResidual()
{
  // Adds - ∫ φ * S(x) dΩ  to the residual
  return -1 * _damage_src[_qp] * _test[_i][_qp];
}

Real
PerturbationSource::computeQpJacobian()
{
  // ∂(φ S)/∂α = 0, since S doesn't depend on α
  return 0.0;
}