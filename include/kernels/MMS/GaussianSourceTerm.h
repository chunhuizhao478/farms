#pragma once
#include "BodyForce.h"

/**
 * Simple Gaussian pulse source for elastic wave propagation with MMS
 */
class GaussianSourceTerm : public BodyForce
{
public:
  static InputParameters validParams();
  GaussianSourceTerm(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  
  // Parameters
  const Real _amplitude;
  const Real _sigma;
  const Real _x0;
  const Real _y0;
  const Real _t0;
  const Real _t_width;
  const Real _omega;
  const unsigned int _component;
};