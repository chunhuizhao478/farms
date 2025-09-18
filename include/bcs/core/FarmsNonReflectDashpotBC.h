#ifndef FarmsNonReflectDashpotBC_H
#define FarmsNonReflectDashpotBC_H

#include "IntegratedBC.h"

/**
 * NonReflecting BC applies a Lysmer damper on a given boundary in the normal
 * and tangential directions
 */
class FarmsNonReflectDashpotBC : public IntegratedBC
{
public:
  static InputParameters validParams();
  FarmsNonReflectDashpotBC(const InputParameters & parameters);

  /**
   * Method for returning parameters that are shared between FarmsNonReflectDashpotBC and
   * FarmsNonReflectDashpotBCAction
   */
  static InputParameters commonParameters();

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Direction in which the Lysmer damper is applied
  unsigned int _component;

  /// Number of displacement variables
  unsigned int _ndisp;

  /// Vector of displacement variables
  std::vector<const VariableValue *> _disp;

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

#endif // FarmsNonReflectDashpotBC_H