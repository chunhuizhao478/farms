//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DiffusedDamageBreakageMaterialMainApp.h"

/**
 * Material used in damage-breakage large deformation formulation for intact porous solids
 * 
 * This material computes effective properties for porous media including:
 * - Effective Biot coefficient
 * - Effective bulk modulus accounting for fluid-solid coupling
 * - Initial elastic parameters for intact (undamaged) material
 * 
 * Inherits damage/breakage framework but represents intact material state
 * Created by Chunhui Zhao, Dec 24th, 2024
 */
class IntactPorousSolidProperties : public Material
{
public:
  static InputParameters validParams();
  IntactPorousSolidProperties(const InputParameters & parameters);

protected:
  /// Compute material properties at quadrature points
  virtual void computeQpProperties() override;

protected:

  Real _lambda_o_value; //lambda_o
  Real _shear_modulus_o_value; //shear_modulus_o
  Real _permeability_solid_o; //permeability_solid_o
  Real _solid_bulk_modulus_s; //solid_bulk_modulus_solid_o
  Real _fluid_bulk_modulus; //fluid_bulk_modulus
  Real _porosity_solid_o; //initial prosoity of solid phase
  Real _initial_viscosity_fluid; // initial viscosity fluid

  MaterialProperty<Real> & _biot_coeff_eff; //_biot_coeff_eff
  MaterialProperty<Real> & _Biot_modulus_eff; //Biot_modulus_eff
  MaterialProperty<Real> & _fluid_solid_coupling; //Fluid_solid_coupling term
  MaterialProperty<Real> & _perm_s; //permeability_solid

};