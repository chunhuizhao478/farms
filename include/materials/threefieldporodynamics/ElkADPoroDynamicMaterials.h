//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADMaterial.h"

/**
 *  ADMaterial used in three field poro dynamics simulations
 */
class ElkADPoroDynamicMaterials : public ADMaterial
{
public:
  static InputParameters validParams();

  ElkADPoroDynamicMaterials(const InputParameters & parameters);

  virtual void computeQpProperties() override;

protected:

  /// solid density value
  ADReal _rhos_val;

  /// fluid density value
  ADReal _rhof_val;

  /// porosity value
  ADReal _porosity_val;

  /// tortosity value
  ADReal _tortosity_val;

  /// viscosity value
  ADReal _viscosity_val;

  /// biot modulus value
  ADReal _bulk_modulus_solid_skeleton_val;

  /// biot coefficient value
  ADReal _bulk_modulus_solid_val;

  /// biot coefficient value
  ADReal _bulk_modulus_fluid_val;

  /// permeability value
  ADReal _permeability_val;

  /// Material property fluid density
  ADMaterialProperty<Real> & _rhof;

  /// Material property density
  ADMaterialProperty<Real> & _rho;

  /// Material property porosity
  ADMaterialProperty<Real> & _porosity;

  /// Material property tortosity
  ADMaterialProperty<Real> & _tortosity;

  /// Material property viscosity
  ADMaterialProperty<Real> & _viscosity;

  /// Material property biot modulus
  ADMaterialProperty<Real> & _biot_modulus;
 
  /// Material property biot coefficients
  ADMaterialProperty<Real> & _biot_coefficient;

  /// Material property permeability
  /// Here the permeability is a scalar
  ADMaterialProperty<Real> & _permeability;
 
};