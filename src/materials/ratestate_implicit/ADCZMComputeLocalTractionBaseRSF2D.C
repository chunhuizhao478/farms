//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADCZMComputeLocalTractionBaseRSF2D.h"
//#include "CZMComputeLocalTractionBase.h"

InputParameters
ADCZMComputeLocalTractionBaseRSF2D::validParams()
{
  InputParameters params = ADCZMComputeLocalTractionBaseRSF2D::validParams();
  return params;
}

ADCZMComputeLocalTractionBaseRSF2D::ADCZMComputeLocalTractionBaseRSF2D(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _interface_traction(
        declareADPropertyByName<RealVectorValue>(_base_name + "interface_traction")),
    _interface_displacement_jump(
        getADMaterialPropertyByName<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _statevar(declareADPropertyByName<Real>("statevar")),
    _sliprate(declareADPropertyByName<Real>("sliprate")),
    _statevarini(getParam<Real>("statevarini")),
    _Vini(getParam<Real>("Vini"))
{
}

void
ADCZMComputeLocalTractionBaseRSF2D::initQpStatefulProperties()
{
  _interface_traction[_qp] = 0;

  _statevar[_qp] = _statevarini;
  _sliprate[_qp] = _Vini;
}

void
ADCZMComputeLocalTractionBaseRSF2D::computeQpProperties()
{
  computeInterfaceTraction();
}
