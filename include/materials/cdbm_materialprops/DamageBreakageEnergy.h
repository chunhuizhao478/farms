//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"


class DamageBreakageEnergy : public Material
{
public:
  static InputParameters validParams();

  DamageBreakageEnergy(const InputParameters & parameters);

  virtual void computeQpProperties() override;

  virtual void initQpStatefulProperties() override;

protected:

  /// Material property initial damage profile
  MaterialProperty<Real> & _Fs; //elastic energy density
  MaterialProperty<Real> & _Fb; //breakage energy density
  MaterialProperty<Real> & _F; //total energy density

  const MaterialProperty<Real> & _lambda;
  const MaterialProperty<Real> & _shear_modulus;
  const MaterialProperty<Real> & _damaged_modulus;
  const MaterialProperty<Real> & _a0;
  const MaterialProperty<Real> & _a1;
  const MaterialProperty<Real> & _a2;
  const MaterialProperty<Real> & _a3;
  const MaterialProperty<Real> & _I1;
  const MaterialProperty<Real> & _I2;
  const MaterialProperty<Real> & _B_breakagevar;

};