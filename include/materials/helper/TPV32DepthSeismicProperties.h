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

/**
 * TPV32DepthSeismicProperties
 *
 * Provides depth-dependent seismic properties (Vs, Vp, rho) following
 * the TPV32 1D piecewise-linear profile. From these, it computes and
 * exposes material properties for shear modulus (mu = rho*Vs^2) and
 * Lamé's first parameter (lambda = rho*Vp^2 - 2*mu).
 *
 * Outputs (MaterialProperty<Real>):
 *   - Vs
 *   - Vp
 *   - rho
 *   - shear_modulus_input   (mu)
 *   - lambda_input          (lambda)
 *
 * Parameters:
 *   - depth_axis: which coordinate represents depth (x/y/z). Default: y.
 *   - depth_offset: value added to the chosen coordinate before lookup.
 *   - clamp_nonpositive_depth: if true (default), depths <= 0 use the 0 m row.
 *   - flip_sign: if true, use depth = -(coord + offset) for negative-down meshes.
 *
 * Output property names (configurable via parameters with these defaults):
 *   - vs_property_name       = "Vs"
 *   - vp_property_name       = "Vp"
 *   - density_property_name  = "density_input"
 *   - shear_modulus_property_name = "shear_modulus_input"
 *   - lambda_property_name   = "lambda_input"
 */
class TPV32DepthSeismicProperties : public Material
{
public:
  static InputParameters validParams();

  TPV32DepthSeismicProperties(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  // Helper: evaluate piecewise linear value at given depth
  Real interp(const std::vector<Real> & xs, const std::vector<Real> & ys, Real x) const;

  // Parameters
  unsigned int _axis; // 0=x,1=y,2=z
  Real _offset;
  bool _clamp_nonpositive;
  bool _flip_sign;

  // Outputs
  MaterialProperty<Real> & _Vs;
  MaterialProperty<Real> & _Vp;
  MaterialProperty<Real> & _rho;
  MaterialProperty<Real> & _mu;
  MaterialProperty<Real> & _lambda;
};
