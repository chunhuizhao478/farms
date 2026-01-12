/*************************************************/
/*           DO NOT MODIFY THIS HEADER           */
/*                                               */
/*                     MASTODON                  */
/*                                               */
/*    (c) 2015 Battelle Energy Alliance, LLC     */
/*            ALL RIGHTS RESERVED                */
/*                                               */
/*   Prepared by Battelle Energy Alliance, LLC   */
/*     With the U. S. Department of Energy       */
/*                                               */
/*     See COPYRIGHT for full restrictions       */
/*************************************************/
#ifndef NonReflectingBCMastodon_H
#define NonReflectingBCMastodon_H

#include "IntegratedBC.h"



/**
 * NonReflecting BC applies a Lysmer damper on a given boundary in the normal
 * and tangential directions
 */
class NonReflectingBCMastodon : public IntegratedBC
{
public:
  static InputParameters validParams();
  NonReflectingBCMastodon(const InputParameters & parameters);

  /**
   * Method for returning parameters that are shared between NonReflectingBCMastodon and
   * NonReflectingBCMastodonAction
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
  std::vector<const VariableValue *> _disp_older;

  /// Density of the soil
  const MaterialProperty<Real> & _density;

  /// P wave speed of the soil
  const Real & _p_wave_speed;

  /// Shear wave speed of the soil
  const Real & _shear_wave_speed;
};

#endif // NonReflectingBCMastodon_H