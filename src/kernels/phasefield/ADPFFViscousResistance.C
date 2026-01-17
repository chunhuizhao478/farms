//* This file is part of the FARMS application
//* Rate-dependent phase field fracture kernel
//* Based on Hofacker & Miehe (2012) - IJNME 93:276-301

#include "ADPFFViscousResistance.h"

registerMooseObject("farmsApp", ADPFFViscousResistance);

InputParameters
ADPFFViscousResistance::validParams()
{
  InputParameters params = ADKernelValue::validParams();
  params.addClassDescription(
      "Viscous resistance term for rate-dependent phase-field fracture. "
      "Implements eta * d_dot = (eta/dt) * (d - d_old) from Hofacker & Miehe (2012). "
      "The weak form is (w, eta * d_dot).");
  params.addParam<MaterialPropertyName>(
      "viscosity", "eta", "The viscosity parameter eta [Pa*s or N*s/m^2]");
  return params;
}

ADPFFViscousResistance::ADPFFViscousResistance(const InputParameters & parameters)
  : ADKernelValue(parameters),
    _eta(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("viscosity"))),
    _d_old(_var.slnOld())
{
}

ADReal
ADPFFViscousResistance::precomputeQpResidual()
{
  // Backward Euler time discretization:
  // d_dot = (d - d_old) / dt
  //
  // Residual contribution to weak form:
  // (w, eta * d_dot) = (w, eta/dt * (d - d_old))
  //
  // For eta = 0, this reduces to the rate-independent case
  // For eta > 0, this provides viscous crack resistance
  return _eta[_qp] * (_u[_qp] - _d_old[_qp]) / _dt;
}
