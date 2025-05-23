#pragma once
#include "CZMComputeLocalTractionTotalBase.h"

/**
 * This material implements a slip-weakening friction law for the method of manufactured solutions
 * based on the 2D dynamic rupture problem with exponential time-dependent displacement.
 */
class MMSSlipWeakeningFriction2dxx : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  MMSSlipWeakeningFriction2dxx(const InputParameters & parameters);
  
protected:
  virtual void computeInterfaceTractionAndDerivatives() override;
  
  // Helper methods for the manufactured solution
  Real spatialFunction(Real x, Real y) const;
  Real temporalFunction(Real t) const;
  Real temporalDerivative(Real t) const;
  Real initialStress(Real x) const;
  Real exactTraction(Real x, Real y, Real t) const;
  
  // Added method to calculate the source term
  Real interfaceSourceTerm(Real x, Real y, Real t, Real sigma_xy, Real sigma_yy) const;
  
  // Coupled variables
  const VariableValue & _nodal_area;
  
  // MMS parameters
  const Real _T2_o;
  const Real _mu_d;
  const Real _Dc;
  const Real _delta;
  const Real _R_m;
  const Real _tbar;
  const Real _tw;
  const Real _Vmin;
  const Real _tau_excess;
  const Real _R_nuc;
  const Real _tau_background;
  const bool _exact_traction;
  
  // Material properties
  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _dstress;
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_old;
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_older;
  const MaterialProperty<Real> & _density;
  
  // Material parameters needed for source term calculation
  const Real _lambda; // First Lamé parameter
  const Real _mu;     // Shear modulus (second Lamé parameter)
  
<<<<<<< HEAD
  // Parameters for the spatial derivatives of the manufactured solution
  const Real _width; // Width parameter for Gaussian spatial function

    // Derived parameters
  static constexpr Real _tau_s = 50e6;    // Static friction strength (Pa)
  static constexpr Real _tau_d = 26.25e6; // Dynamic friction strength (Pa)
=======
  // Derived parameters
  static constexpr Real _tau_s = 50e6;    // Static friction strength (Pa)
  static constexpr Real _tau_d = 26.25e6; // Dynamic friction strength (Pa)
  
  // Parameters for the spatial derivatives of the manufactured solution
  const Real _width; // Width parameter for Gaussian spatial function
>>>>>>> 6188d19945ce13b4debb058d2e709731e3f73bb8
};