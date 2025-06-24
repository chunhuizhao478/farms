/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Gaussian Distribution
Chunhui Zhao
*/

#pragma once

#include "Function.h"

class InitialStressStrainCDBMv2 : public Function
{
public:
  InitialStressStrainCDBMv2(const InputParameters & parameters);

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
  Real _peak_shear_value; //initial shear stress perturbation peak value
  Real _nucl_center_x; //nucleation center x coordinate
  Real _nucl_center_z; //nucleation center z coordinate
  Real _nucl_size; //nucleation size
  Real _elem_size; //element size for the simulation, used for determining the nucleation zone
  Real _cutoff_distance; //cutoff distance for the depth varying stress
  bool _get_initial_stress; //flag to get initial stress
  bool _get_initial_strain; //flag to get initial strain
  bool _get_shear_overstress; //flag to get initial shear overstress
  bool _get_fluid_pressure; //flag to get fluid pressure

};