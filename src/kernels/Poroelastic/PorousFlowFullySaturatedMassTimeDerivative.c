#include "PorousFlowFullySaturatedMassTimeDerivative.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableFE.h"
#include "libmesh/quadrature.h"

registerMooseObject("farmsApp", PorousFlowFullySaturatedMassTimeDerivative);

InputParameters
PorousFlowFullySaturatedMassTimeDerivative::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addClassDescription("The time derivative operator with the weak form of $(\\psi_i, "
                             "\\frac{\\partial u_h}{\\partial t})$.");
  params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");
  return params;
}

PorousFlowFullySaturatedMassTimeDerivative::PorousFlowFullySaturatedMassTimeDerivative(const InputParameters & parameters)
  : TimeKernel(parameters), _lumping(getParam<bool>("lumping")),
     _biot_modulus(getMaterialProperty<Real>("PorousFlow_constant_biot_modulus_qp")),
{
}

Real
PorousFlowFullySaturatedMassTimeDerivative::computeQpResidual()
{
  return _test[_i][_qp] * _u_dot[_qp]/_biot_modulus[_qp];
}

Real
PorousFlowFullySaturatedMassTimeDerivative::computeQpJacobian()
{
  return _test[_i][_qp] * _phi[_j][_qp] * _du_dot_du[_qp]/_biot_modulus[_qp];
}

void
PorousFlowFullySaturatedMassTimeDerivative::computeJacobian()
{
  if (_lumping)
  {
    prepareMatrixTag(_assembly, _var.number(), _var.number());

    precalculateJacobian();
    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < _phi.size(); _j++)
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
          _local_ke(_i, _i) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();

    accumulateTaggedLocalMatrix();
  }
  else
    TimeKernel::computeJacobian();
}