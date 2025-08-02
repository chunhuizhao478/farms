#pragma once

#include "ADIntegratedBC.h"

/**
 * NonReflecting BC applies a Lysmer damper on a given boundary in the normal
 * and tangential directions using automatic differentiation
 */
class ADFarmsNonReflectDashpotBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  ADFarmsNonReflectDashpotBC(const InputParameters & parameters);

  /**
   * Method for returning parameters that are shared between ADFarmsNonReflectDashpotBC and
   * ADFarmsNonReflectDashpotBCAction
   */
  static InputParameters commonParameters();

protected:
  virtual ADReal computeQpResidual() override;

  /// Direction in which the Lysmer damper is applied
  const unsigned int _component;

  /// Number of displacement variables
  const unsigned int _ndisp;

  /// Vector of displacement variables
  std::vector<const ADVariableValue *> _disp;

  /// Unsigned integers representing the displacement variables
  std::vector<unsigned int> _disp_var;

  /// Vector of old displacement variables
  std::vector<const VariableValue *> _disp_old;

  /// Vector of old velocity variables
  std::vector<const VariableValue *> _vel_old;

  /// Vector of old acceleration variables
  std::vector<const VariableValue *> _accel_old;

  /// _beta Parameter for Newmark time integration scheme
  const Real _beta;

  /// _gamma Parameter for Newmark time integration scheme
  const Real _gamma;

  /// _alpha Parameter for HHT time integration scheme
  const Real _alpha;

  /// Density of the soil
  const Real & _density;

  /// P wave speed of the soil
  const Real & _p_wave_speed;

  /// Shear wave speed of the soil
  const Real & _shear_wave_speed;
};