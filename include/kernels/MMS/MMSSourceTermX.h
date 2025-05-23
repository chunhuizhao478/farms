#pragma once

#include "BodyForce.h"

/**
 * MMSSourceTermX implements the source term for Method of Manufactured Solutions 
 * in x-direction for a continuous domain (without slip weakening)
 */
class MMSSourceTermX : public BodyForce
{
public:
  static InputParameters validParams();
  MMSSourceTermX(const InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual() override;
  
  // Helper functions for manufactured solution
  Real spatialFunction(Real x, Real y) const;
  Real temporalFunction(Real t) const;
  Real temporalDerivative(Real t) const;
  Real temporalSecondDerivative(Real t) const;
  Real sourceTerm(Real x, Real y, Real t) const;
  
  // Parameters for manufactured solution
  const Real _delta;              // Displacement amplitude
  const Real _R_m;                // Characteristic length (width = _R_m/2)
  const Real _tbar;               // Rupture time
  const Real _tw;                 // Event duration parameter
  const Real _Vmin;               // Minimum velocity
  
  // Material properties
  const Real _density;            // Density
  const Real _lambda;             // First Lamé parameter
  const Real _mu;                 // Shear modulus (second Lamé parameter)
  
  // Background stress
  const Real _tau_background;     // Background stress
};