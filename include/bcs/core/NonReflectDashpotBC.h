//Implement Lysmer damper ( same as NonreflectingBC | MASTODON ) but in explicit scheme ( follow DashpotBC | MOOSE )
//Created by Chunhui Zhao, Jan 18, 2023

#pragma once

#include "IntegratedBC.h"

class NonReflectDashpotBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  NonReflectDashpotBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

private:
  /// Component of the velocity vector
  unsigned int _component;

  unsigned int _disp_x_var;
  unsigned int _disp_y_var;

  const VariableValue & _disp_x_dot;
  const VariableValue & _disp_y_dot;

  const MaterialProperty<Real> & _density;

  const Real & _p_wave_speed;
  const Real & _shear_wave_speed;
};