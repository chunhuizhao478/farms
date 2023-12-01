/*
Define pore pressure distribution given analytical solution
Reference: John W. RUDNICKI : FLUID MASS SOURCES AND POINT FORCES IN LINEAR ELASTIC DIFFUSIVE SOLIDS

Chunhui Zhao
*/

#pragma once

#include "Function.h"

class InitialStressXYPressureNetwork : public Function
{
public:
  InitialStressXYPressureNetwork(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  Real _flux_q;            //flux
  Real _density_rho_0;     //fluid density
  Real _permeability_k;    //permeability
  Real _viscosity_eta;     //viscosity
  Real _biotcoeff_alpha;   //biot coefficient
  Real _undrained_nu_u;    //undrained poisson's ratio
  Real _shear_modulus_mu;  //shear modulus
  Real _drained_nu;        //drained poisson's ratio

};