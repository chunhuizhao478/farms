// In AdaptiveTimeStepCalculator.h
#pragma once

#include "TimeStepper.h"
#include "PostprocessorInterface.h"

class AdaptiveTimeStepCalculator : public TimeStepper,
                                  public PostprocessorInterface
{
public:
  static InputParameters validParams();
  AdaptiveTimeStepCalculator(const InputParameters & parameters);

protected:
  virtual void init() override;
  virtual Real computeInitialDT() override;
  virtual Real computeDT() override;
  
private:
  // Parameters from input file
  const Real _cp;
  const Real _a_prem;
  const Real _shear_modulus;
  const Real _permeability;
  const Real _biot_modulus;
  
  const Real _a_o;
  const Real _b_o;
  const Real _L;
  const Real _normal_stress;
  const Real _zeta_max;

  const PostprocessorValue & _dx_min;
  const PostprocessorValue & _max_slip_rate;
  
  // Calculated parameters
  Real _dt_exp;
  Real _dt_diff;
  Real _zeta;
};