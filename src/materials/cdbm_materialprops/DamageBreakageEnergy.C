//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DamageBreakageEnergy.h"

registerMooseObject("farmsApp", DamageBreakageEnergy);

InputParameters
DamageBreakageEnergy::validParams()
{
  InputParameters params = Material::validParams();
  return params;
}

DamageBreakageEnergy::DamageBreakageEnergy(const InputParameters & parameters)
  : Material(parameters),
  //declare properties
  _Fs(declareProperty<Real>("elastic_energy_density")),
  _Fb(declareProperty<Real>("breakage_energy_density")),
  _F(declareProperty<Real>("total_energy_density")),
  _lambda(getMaterialProperty<Real>("lambda_const")),
  _shear_modulus(getMaterialProperty<Real>("shear_modulus")),
  _damaged_modulus(getMaterialProperty<Real>("damaged_modulus")),
  _a0(getMaterialProperty<Real>("a0")),
  _a1(getMaterialProperty<Real>("a1")),
  _a2(getMaterialProperty<Real>("a2")),
  _a3(getMaterialProperty<Real>("a3")),
  _I1(getMaterialProperty<Real>("first_elastic_strain_invariant")),
  _I2(getMaterialProperty<Real>("second_elastic_strain_invariant")),
  _B_breakagevar(getMaterialProperty<Real>("B_damagedvar"))
{
}

void
DamageBreakageEnergy::initQpStatefulProperties()
{
}

void
DamageBreakageEnergy::computeQpProperties()
{

  //compute elastic energy density
  _Fs[_qp] = 0.5 * _lambda[_qp] * _I1[_qp] * _I1[_qp] + _shear_modulus[_qp] * _I2[_qp] - _damaged_modulus[_qp] * _I1[_qp] * sqrt(_I2[_qp]);

  //compute breakage energy density
  _Fb[_qp] = _a0[_qp] * _I2[_qp] + _a1[_qp] * _I1[_qp] * sqrt(_I2[_qp]) + _a2[_qp] * _I2[_qp] * _I2[_qp] + _a3[_qp] * _I1[_qp] * _I1[_qp] * _I1[_qp] / sqrt(_I2[_qp]);

  //compute total energy density
  _F[_qp] = _Fs[_qp] + _Fb[_qp];

}