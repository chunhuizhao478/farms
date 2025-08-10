//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkADVarDivTest.h"
#include "Assembly.h"

registerMooseObject("farmsApp", ElkADVarDivTest);

InputParameters
ElkADVarDivTest::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The gradient operator optionally scaled by a constant scalar "
                             "coefficient. Weak form: $(\\nabla \\cdot \\vec{\\psi_i}, k v)$.");
  params.addRequiredCoupledVar("coupled_scalar_variable", "The scalar field");
  params.addRequiredCoupledVar("w_x", "The darcy velocity along x direction");
  params.addRequiredCoupledVar("w_y", "The darcy velocity along y direction");
  return params;
}

ElkADVarDivTest::ElkADVarDivTest(const InputParameters & parameters)
  : ADKernel(parameters),
    _p(adCoupledValue("coupled_scalar_variable")),
    _wx(getVar("w_x", 0)),
    _wy(getVar("w_y", 0)),
    _grad_wx_phi(_wx->gradPhi()),
    _grad_wy_phi(_wy->gradPhi())
{
}

ADReal
ElkADVarDivTest::computeQpResidual()
{
  
  ADReal div_phi = 0;
  div_phi = _grad_wx_phi[_i][_qp](0) + _grad_wy_phi[_i][_qp](1);

  return - div_phi * _p[_qp];
}
