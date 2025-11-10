/*
Function to compute initial stress and strain based on SCEC TPV26 Benchmark.
Created By Chunhui Zhao, Jul 11th, 2025
*/

#pragma once

#include "Function.h"

class InitialStressStrainTPV26 : public Function
{
public:
  InitialStressStrainTPV26(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  Real _i; //index
  Real _j; //index
  Real _lambda_o; //initial lambda parameter for the CDBM model
  Real _shear_modulus_o; //initial shear modulus parameter for the CDBM model
  Real _fluid_density; //fluid density in kg/m^3
  Real _rock_density; //rock density in kg/m^3
  Real _gravity; //gravity in m/s^2
  Real _bxx; //coefficient for sigmaxx
  Real _byy; //coefficient for sigmayy
  Real _bxy; //coefficient for sigmaxy
  bool _get_initial_stress; //flag to get initial stress
  bool _get_initial_strain; //flag to get initial strain
  bool _get_shear_overstress; //flag to get initial shear overstress
  bool _get_fluid_pressure; //flag to get fluid pressure
  bool _use_tapering; //flag to use tapering in the stress calculation to reduce deviatroic stress components start at a certain depth, default is false
  Real _tapering_depth_A; //depth at which tapering starts to be applied
  Real _tapering_depth_B; //depth at which tapering stops to be applied
  bool _use_overpressure; //flag to use overpressure in the stress calculation
  Real _overpressure_depth_A; //depth at which overpressure starts to transition from hydrostatic to lithostatic
  Real _overpressure_depth_B; //depth at which overpressure stops to transition from hydrostatic to lithostatic
  bool _overpressure_loweffective; //flag to use low effective stress overpressure (quadratic transition, lambda_pp scaling)
  Real _lambda_pp; //pore pressure ratio for low effective stress overpressure (Pf = lambda_pp * rho * g * z below depth B)

};