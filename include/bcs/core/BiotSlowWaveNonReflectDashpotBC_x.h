//Implement Lysmer damper ( same as NonreflectingBC | MASTODON ) but in explicit scheme ( follow DashpotBC | MOOSE )
//Created by Amr Ibrahim, April 18, 2024

#pragma once

#include "IntegratedBC.h"

class BiotSlowWaveNonReflectDashpotBC_x : public IntegratedBC
{
public:
  static InputParameters validParams();

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  BiotSlowWaveNonReflectDashpotBC_x(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  /// Component of the velocity vector
  unsigned int _component;

  const VariableValue & _fluid_vel_y;

  unsigned int _fluid_vel_y_id;

  const VariableValue & _fluid_vel_z;

  unsigned int _fluid_vel_z_id;

  const MaterialProperty<Real> & _rho_f;

  const Real & _biot_p_wave_speed;
};