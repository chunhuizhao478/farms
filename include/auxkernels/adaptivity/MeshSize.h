//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

/**
 * AuxKernel that computes current elemental size 
 */
class MeshSize : public AuxKernel
{
public:
  static InputParameters validParams();

  MeshSize(const InputParameters & parameters);

protected:

  virtual Real computeValue() override;

  //Energy release rate as constant material property
  Real _energy_release_rate; 

  //Young's Modulus
  Real _youngs_modulus;

  //Weibull distribution perturbation parameters
  Real _seed;
  Real _shape;
  Real _scale;
  Real _min;
  Real _max;

};