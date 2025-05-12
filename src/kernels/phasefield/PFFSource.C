//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

//Note:
/**
 * The KernelValue class is responsible for calculating the residuals in form:
 *
 *  JxW[_qp] * _value[_qp] * _test[_i][_qp]
 *
 */

#include "PFFSource.h"

registerMooseObject("farmsApp", PFFSource);

InputParameters
PFFSource::validParams()
{
  InputParameters params = KernelValue::validParams();
  params.addClassDescription("The source term in the phase-field evolution equation. The weak form "
                             "is $(w, \\dfrac{\\partial \\psi}{\\partial d})$.");
  params.addParam<MaterialPropertyName>("free_energy_first_derivative",
                                        "dpsi_dd",
                                        "The first derivative of Helmholtz free energy");
  params.addParam<MaterialPropertyName>("free_energy_second_derivative",
                                        "d2psi_dd2",
                                        "The second derivative of Helmholtz free energy");
  return params;
}

PFFSource::PFFSource(const InputParameters & parameters)
  : KernelValue(parameters),
    DerivativeMaterialPropertyNameInterface(),
    _dpsi_dd(getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("free_energy_first_derivative"))),
    _d2psi_dd2(getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("free_energy_second_derivative")))
{
}

Real
PFFSource::precomputeQpResidual()
{
  return _dpsi_dd[_qp];
}

Real
PFFSource::precomputeQpJacobian()
{
  return _d2psi_dd2[_qp];
}