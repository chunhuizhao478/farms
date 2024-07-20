//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADComputeEigenstrainBase.h"
#include "RankFourTensor.h"

/**
 * ComputeEigenstrain computes an Eigenstrain
 */
class ADComputeEigenstrainFromSolution : public ADComputeEigenstrainBase
{
public:
  static InputParameters validParams();

  ADComputeEigenstrainFromSolution(const InputParameters & parameters);

protected:
  virtual void computeQpEigenstrain() override;

  /// base_name for elasticity tensor to use to convert strain to strain
  const std::string _base_name;

  ///Stores the total eigenstrain in the previous step
  const MaterialProperty<RankTwoTensor> & _eigenstrain_old;

  /// Whether the user has supplied AuxVariables representing the initial strain
  const bool _ini_aux_provided;

  /// initial strain components
  std::vector<const Function *> _initial_strain_fcn;

  /// AuxVariables defining the initial strain
  const std::vector<const VariableValue *> _ini_aux;
};
