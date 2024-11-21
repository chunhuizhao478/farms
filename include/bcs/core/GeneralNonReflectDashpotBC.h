//Implement Lysmer damper ( same as NonreflectingBC | MASTODON ) but in explicit scheme ( follow DashpotBC | MOOSE )
//Created by Amr Ibrahim, April 18, 2024

#pragma once

#include "IntegratedBC.h"

class GeneralNonReflectDashpotBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  GeneralNonReflectDashpotBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  /// Component of the velocity vector
  unsigned int _component;

  /// number of displacements
  const unsigned int _ndisp;

  /// Variable numbers of the displacement variables
  const std::vector<unsigned int> _disp_var;


  const VectorVariableValue & _disp_dot;
  const VariableValue & _dvel_dot;
  const VariableValue & _du_dot_du;

  const MaterialProperty<Real> & _density;

  const Real & _p_wave_speed;
  const Real & _shear_wave_speed;
};
