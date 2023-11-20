#include "PorousFlowFullySaturatedVolumetricExpansion.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableFE.h"
#include "libmesh/quadrature.h"

registerMooseObject("farmsApp", PorousFlowFullySaturatedVolumetricExpansion);

InputParameters
PorousFlowFullySaturatedVolumetricExpansion::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addClassDescription("The time derivative operator with the weak form of $(\\psi_i, "
                             "\\frac{\\partial u_h}{\\partial t})$.");
  params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");

  params.addRangeCheckedParam<Real>(
      "biot_coefficient", 1.0, "biot_coefficient>=0 & biot_coefficient<=1", "Biot coefficient");

  return params;
}

PorousFlowFullySaturatedVolumetricExpansion::PorousFlowFullySaturatedVolumetricExpansion(const InputParameters & parameters)
  : TimeKernel(parameters), _lumping(getParam<bool>("lumping"))
{
}

Real
PorousFlowFullySaturatedVolumetricExpansion::computeQpResidual()
{
  return _test[_i][_qp]* _biot_coefficient* _u_dot[_qp]/3;
}

Real
PorousFlowFullySaturatedVolumetricExpansion::computeQpJacobian()
{
  return _test[_i][_qp] * _biot_coefficient* _phi[_j][_qp] * _du_dot_du[_qp]/3;
}

void
PorousFlowFullySaturatedVolumetricExpansion::computeJacobian()
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