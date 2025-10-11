#pragma once

#include "Material.h"
#include "MooseEnum.h"

/**
 * Material that computes the spatially varying attenuation coefficient
 * used by the Sochacki sponge boundary. The coefficient is built from a
 * normalized distance into the sponge region and applies one of the
 * damping profiles (linear, polynomial, exponential, Gaussian) described
 * by Sochacki et al. (1987).
 */
class SochackiSpongeMaterial : public Material
{
public:
  static InputParameters validParams();
  SochackiSpongeMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

private:
  /// Helper that returns normalized distance (0 at interface, 1 at outer edge)
  Real computeNormalizedDistance(const Point & p) const;

  MaterialProperty<Real> & _sponge_coeff;

  const Real _inner_xmin;
  const Real _inner_xmax;
  const Real _inner_ymin;
  const Real _inner_ymax;

  const Real _outer_xmin;
  const Real _outer_xmax;
  const Real _outer_ymin;
  const Real _outer_ymax;

  const Real _s_max;
  const Real _exponent_power;
  const Real _exp_rate;
  const Real _gaussian_rate;

  const MooseEnum _profile;

  const Real _left_thickness;
  const Real _right_thickness;
  const Real _bottom_thickness;
  const Real _top_thickness;

  const Real _exp_norm;
  const Real _gaussian_norm;
};
