#include "InitialShearTractionKernel.h"

registerMooseObject("farmsApp", InitialShearTractionKernel);

InputParameters
InitialShearTractionKernel::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addRequiredCoupledVar("ini_shear_sts_perturb","initial shear stress perturbation spatial distribution");
  return params;
}

InitialShearTractionKernel::InitialShearTractionKernel(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _ini_shear_sts_perturb(coupledValue("ini_shear_sts_perturb"))
{
}

Real
InitialShearTractionKernel::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::Element:
      r = _test[_i][_qp] * _ini_shear_sts_perturb[_qp];
      break;

    case Moose::Neighbor:
      r = -_test_neighbor[_i][_qp] * _ini_shear_sts_perturb[_qp];
      break;
  }

  return r;
}