//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "PFFDiffusion.h"

registerMooseObject("farmsApp", PFFDiffusion);

InputParameters
PFFDiffusion::validParams()
{
  InputParameters params = Kernel::validParams();
  params += BaseNameInterface::validParams();
  params.addClassDescription("The diffusion term in the phase-field evolution equation. The weak "
                             "form is $(\\grad w, \\dfrac{2\\Gc l}{c_0} \\grad d)$.");

  params.addParam<MaterialPropertyName>(
      "fracture_toughness", "Gc", "The fracture toughness $\\Gc$");
  params.addParam<MaterialPropertyName>(
      "normalization_constant", "c0", "The normalization constant $c_0$");
  params.addParam<MaterialPropertyName>(
      "regularization_length", "l", "The phase-field regularization length");
  return params;
}

PFFDiffusion::PFFDiffusion(const InputParameters & parameters)
  : Kernel(parameters),
    BaseNameInterface(parameters),
    _Gc(getMaterialProperty<Real>(prependBaseName("fracture_toughness", true))),
    _c0(getMaterialProperty<Real>(prependBaseName("normalization_constant", true))),
    _l(getMaterialProperty<Real>(prependBaseName("regularization_length", true)))
{
}

Real
PFFDiffusion::computeQpResidual()
{
  Real value = _grad_test[_i][_qp] * _grad_u[_qp];

  return 2 * _Gc[_qp] * _l[_qp] / _c0[_qp] * value;
}

Real
PFFDiffusion::computeQpJacobian()
{
  Real coef = 2 * _Gc[_qp] * _l[_qp] / _c0[_qp];

  //same as diffusion.C
  return coef * _grad_test[_i][_qp] * _grad_phi[_j][_qp];
}
