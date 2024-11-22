// SwitchStateFunction.h
#pragma once

#include "Function.h"
#include "MaterialProperty.h"

class SwitchStateFunction : public Function
{
public:
  static InputParameters validParams();
  SwitchStateFunction(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) const override;

private:
  const PostprocessorValue & _vel_sol;
  const Real _vel_threshold_quasi_to_dyn;
  const Real _vel_threshold_dyn_to_quasi;
  const PostprocessorValue & _stateflag_old;
  const PostprocessorValue & _stateflag_older;
  const PostprocessorValue & _stateflag_oldolder;
};