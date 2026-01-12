/*
initial Cohesion function for CDBMv2 benchmark, defined by a piecewise linear function based on depth.
Created By Chunhui Zhao, Jun 24th, 2025
*/

#pragma once

#include "Function.h"

class InitialCohesionCDBMv2 : public Function
{
public:
  InitialCohesionCDBMv2(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  Real _depth;        // Depth at which cohesion starts to decrease (in meters)
  Real _slope;        // Slope of the linear decrease in cohesion with depth (in MPa per meter)
  Real _min_cohesion; // Minimum cohesion value at the maximum depth (in MPa)

};