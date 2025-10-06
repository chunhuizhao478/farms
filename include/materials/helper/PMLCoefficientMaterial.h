//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*
 * Material that computes PML damping coefficient
 * The damping increases gradually from the physical domain into the PML layer
 */

#pragma once

#include "Material.h"

/**
 * PMLCoefficientMaterial computes the spatially-varying damping coefficient
 * for the Perfectly Matched Layer (PML).
 *
 * The damping profile follows: d(r) = d_max * (r/L)^n
 * where:
 *   r = distance from PML inner boundary
 *   L = PML thickness
 *   n = damping exponent (typically 2-4)
 *   d_max = maximum damping coefficient
 */
class PMLCoefficientMaterial : public Material
{
public:
  static InputParameters validParams();

  PMLCoefficientMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

private:
  /// PML damping coefficient (output property)
  MaterialProperty<Real> & _pml_damping_coeff;

  /// Inner boundary of PML (physical domain edge)
  const Real _pml_xmin;
  const Real _pml_xmax;
  const Real _pml_ymin;
  const Real _pml_ymax;

  /// PML thickness
  const Real _pml_thickness;

  /// Maximum damping coefficient
  const Real _d_max;

  /// Damping profile exponent
  const Real _exponent;

  /// Reference wave speed (for scaling d_max)
  const Real _ref_wave_speed;

  /**
   * Compute distance into PML region
   * Returns 0 in physical domain, increases in PML
   */
  Real computePMLDistance(const Point & p) const;
};
