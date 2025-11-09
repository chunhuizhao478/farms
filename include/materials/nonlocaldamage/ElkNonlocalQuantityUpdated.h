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
#include "ElkRadialAverageUpdated.h"

/**
 * ElkNonlocalQuantityUpdated retrieves a nonlocal averaged material property
 * from ElkRadialAverageUpdated user object.
 *
 * This is a generic class that can handle any scalar material property that
 * needs nonlocal averaging (e.g., strain invariant ratio, strain rate, etc.)
 *
 * Usage:
 *   For strain invariant ratio (backward compatible):
 *     type = ElkNonlocalQuantityUpdated
 *     average_UO = eqstrain_averaging_block1
 *     nonlocal_property_name = 'eqstrain_nonlocal'
 *     initial_property_name = 'eqstrain_nonlocal_initial'
 *
 *   For strain rate:
 *     type = ElkNonlocalQuantityUpdated
 *     average_UO = strainrate_averaging_block1
 *     nonlocal_property_name = 'strain_rate_nonlocal'
 *     initial_property_name = 'strain_rate_nonlocal_initial'
 */
class ElkNonlocalQuantityUpdated : public Material
{
public:
  static InputParameters validParams();

  ElkNonlocalQuantityUpdated(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// The nonlocal averaged material property (declared by this material)
  MaterialProperty<Real> & _prop_nonlocal;

  /// Old value of nonlocal property
  const MaterialProperty<Real> & _prop_nonlocal_old;

  /// Reference to the averaging user object
  const std::map<dof_id_type, std::vector<Real>> & _average;

  /// Cache for element lookup optimization
  const Elem * _prev_elem;
  std::map<dof_id_type, std::vector<Real>>::const_iterator _average_iterator;

  /// Initial value property
  const MaterialProperty<Real> & _prop_initial;

  /// Time step
  int & _step;
};
