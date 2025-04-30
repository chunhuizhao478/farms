/*
Define Function for Experiement Loading Setup
Created by Chunhui Zhao, Oct 26th, 2024
*/

#pragma once

#include "Function.h"

class ElkPulseLoadExperiment : public Function
{
public:
  ElkPulseLoadExperiment(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  Real _shape_param_alpha;
  Real _shape_param_beta;
  Real _rise_time;
  Real _single_pulse_duration;
  Real _convert_efficiency;
  Real _EM;
  Real _gap;
  Real _fitting_param_alpha;
  std::vector<Real> _discharge_center;
  int  _number_of_pulses;

  Real _peak_pressure;

  //minimum applied pressure applied on the boundary
  //to mimic the effect of water pressure
  Real _minimum_applied_pressure;
  bool _use_minimum_applied_pressure;

};