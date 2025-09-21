/*
Define Function for Experiement Loading Setup
Created by Chunhui Zhao, Sep 21th, 2025
Pmax is determined following Wu 2022 paper
*/

#pragma once

#include "Function.h"

class ElkPulseLoadExperimentWu2022Paper : public Function
{
public:
  ElkPulseLoadExperimentWu2022Paper(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  Real _shape_param_alpha;
  Real _shape_param_beta;
  Real _rise_time;
  Real _single_pulse_duration;
  std::vector<Real> _discharge_center;
  std::vector<Real> _pmax_coefficients;
  int  _number_of_pulses;

  //minimum applied pressure applied on the boundary
  //to mimic the effect of water pressure
  Real _minimum_applied_pressure;
  bool _use_minimum_applied_pressure;

};