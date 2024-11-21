//Implement Lysmer damper ( same as NonreflectingBC | MASTODON ) but in explicit scheme ( follow DashpotBC | MOOSE )
//Created by Amr Ibrahim, April 18, 2024

#pragma once

#include "IntegratedBC.h"

class GeneralNonReflectDashpotBC_y : public IntegratedBC
{
public:
  static InputParameters validParams();

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  GeneralNonReflectDashpotBC_y(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  /// Component of the velocity vector
  unsigned int _component;

  const VariableValue & _disp_x_dot;
  const VariableValue & _d_disp_x_dot;

  unsigned int _disp_x_id;

  const VariableValue & _disp_z_dot;
  const VariableValue & _d_disp_z_dot;

  unsigned int _disp_z_id;

  const VariableValue & _u_dot;
  const VariableValue & _du_dot_du;


  const MaterialProperty<Real> & _density;

  const Real & _p_wave_speed;
  const Real & _shear_wave_speed;
};