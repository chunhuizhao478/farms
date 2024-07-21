//Implement Lysmer damper ( same as NonreflectingBC | MASTODON ) but in explicit scheme ( follow DashpotBC | MOOSE )
//Generalization to 3D
//Created by Chunhui Zhao, May 1, 2023

#pragma once

#include "ADIntegratedBC.h"

class ADNonReflectDashpotBC3d : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  ADNonReflectDashpotBC3d(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

private:
  /// Component of the velocity vector
  unsigned int _component;

  unsigned int _disp_x_var;
  unsigned int _disp_y_var;
  unsigned int _disp_z_var;

  const ADVariableValue & _disp_x_dot;
  const ADVariableValue & _disp_y_dot;
  const ADVariableValue & _disp_z_dot;

  const ADMaterialProperty<Real> & _density;

  const Real & _p_wave_speed;
  const Real & _shear_wave_speed;
};