//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SochackiSpongeDamping.h"

registerMooseObject("farmsApp", SochackiSpongeDamping);

InputParameters
SochackiSpongeDamping::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Implements the damping term 2*A(x,y)*rho*du/dt used in the "
                             "Sochacki sponge absorbing boundary.");
  params.addParam<Real>("damping_scale",
                        2.0,
                        "Multiplicative scale applied to the attenuation coefficient. "
                        "The canonical value for the Sochacki sponge is 2.0.");
  params.addRequiredParam<MaterialPropertyName>(
      "density", "Name of the density material property used in the momentum balance.");
  params.addRequiredParam<MaterialPropertyName>(
      "sochacki_damping", "Attenuation coefficient supplied by SochackiSpongeMaterial.");
  return params;
}

SochackiSpongeDamping::SochackiSpongeDamping(const InputParameters & parameters)
  : Kernel(parameters),
    _u_dot(dot()),
    _density(getMaterialProperty<Real>(getParam<MaterialPropertyName>("density"))),
    _sponge_coeff(getMaterialProperty<Real>(getParam<MaterialPropertyName>("sochacki_damping"))),
    _damping_scale(getParam<Real>("damping_scale"))
{
}

Real
SochackiSpongeDamping::computeQpResidual()
{
  return _test[_i][_qp] * _damping_scale * _sponge_coeff[_qp] * _density[_qp] * _u_dot[_qp];
}

Real
SochackiSpongeDamping::computeQpJacobian()
{
  // Damping term depends on du/dt; for explicit central difference the Jacobian is zero.
  return 0.0;
}
