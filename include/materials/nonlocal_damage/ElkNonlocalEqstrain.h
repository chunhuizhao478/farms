//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"
#include "Material.h"
#include "ElkRadialAverage.h"

// Forward declarations

/*
Declare a nonlocal eqstrain material property ""
update it using "ElkRadialAverage"
Then the nonlocal eqstrain material property will be in the initial values of "eqstrain_nonlocal"
in the "ElkComputeSmearedCrackingStressModifiedMazars"

Created By Chunhui Zhao, Oct 19th, 2024
*/

class ElkNonlocalEqstrain : public Material
{
public:
  static InputParameters validParams();

  ElkNonlocalEqstrain(const InputParameters & parameters);

  virtual void initQpStatefulProperties() override;

  virtual void computeQpProperties() override;

protected:

  /// @brief define nonlocal eqstrain
  MaterialProperty<Real> & _eqstrain_nonlocal;
  const MaterialProperty<Real> & _eqstrain_nonlocal_old;

  // Averaged Material
  // std::string _avg_material_name;
  const ElkRadialAverage::Result & _average;
  ElkRadialAverage::Result::const_iterator _average_eqstrain_nonlocal;

  // Pointer to last element for comparison for speed
  const Elem * _prev_elem;

};