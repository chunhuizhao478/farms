//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"
#include "RankTwoTensor.h"

// Forward declarations

class ADElasticEnergyAux : public AuxKernel
{
public:
  static InputParameters validParams();

  ADElasticEnergyAux(const InputParameters & parameters);
  virtual ~ADElasticEnergyAux() {}

protected:
  virtual Real computeValue();

  /// Base name of the material system used to calculate the elastic energy
  const std::string _base_name;

  /// The stress tensor
  const ADMaterialProperty<RankTwoTensor> & _stress;
  const ADMaterialProperty<RankTwoTensor> & _elastic_strain;
};