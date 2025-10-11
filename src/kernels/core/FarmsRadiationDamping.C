//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FarmsRadiationDamping.h"

registerMooseObject("farmsApp", FarmsRadiationDamping);

InputParameters
FarmsRadiationDamping::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addClassDescription("Applies domain-wide radiation damping proportional to velocity.");
  params.addRangeCheckedParam<Real>(
      "eta_constant",
      0.0,
      "eta_constant>=0.0",
      "Constant damping coefficient multiplying the velocity term.");
  params.addParam<MaterialPropertyName>(
      "eta_property",
      "Optional material property providing spatially varying damping coefficient.");
  params.addParam<MaterialPropertyName>(
      "flag_property",
      "Optional material property that scales the damping contribution, e.g. a 0/1 switch.");

  return params;
}

FarmsRadiationDamping::FarmsRadiationDamping(const InputParameters & parameters)
  : TimeDerivative(parameters),
    _eta_constant(getParam<Real>("eta_constant")),
    _eta_property(parameters.isParamSetByUser("eta_property")
                      ? &getMaterialProperty<Real>(getParam<MaterialPropertyName>("eta_property"))
                      : nullptr),
    _flag_property(parameters.isParamSetByUser("flag_property")
                       ? &getMaterialProperty<Real>(getParam<MaterialPropertyName>("flag_property"))
                       : nullptr)
{
  if (!_eta_property && !parameters.isParamSetByUser("eta_constant"))
    paramError("eta_constant",
               "Either 'eta_constant' must be provided or 'eta_property' must be specified.");
}

Real
FarmsRadiationDamping::dampingCoefficient() const
{
  Real eta = _eta_property ? (*_eta_property)[_qp] : _eta_constant;
  if (_flag_property)
    eta *= (*_flag_property)[_qp];

  return eta;
}

Real
FarmsRadiationDamping::computeQpResidual()
{
  return dampingCoefficient() * TimeDerivative::computeQpResidual();
}

Real
FarmsRadiationDamping::computeQpJacobian()
{
  return dampingCoefficient() * TimeDerivative::computeQpJacobian();
}

Real
FarmsRadiationDamping::computeQpOffDiagJacobian(unsigned int jvar)
{
  return dampingCoefficient() * TimeDerivative::computeQpOffDiagJacobian(jvar);
}
