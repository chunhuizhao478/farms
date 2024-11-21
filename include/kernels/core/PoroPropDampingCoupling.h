#pragma once

#include "Kernel.h"

class PoroPropDampingCoupling : public Kernel
{
public:
  static InputParameters validParams();

  PoroPropDampingCoupling(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  /// An integer corresponding to the direction this kernel acts in
  unsigned int _component;

  /// Biot coefficient
  const MaterialProperty<Real> & _coefficient;

  const VariableValue & _porepressure;

  const VariableValue & _porepressure_old;

  unsigned int _porepressure_var_num;

  Real _q;

};