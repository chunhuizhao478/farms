#pragma once
#include "Kernel.h"

/**
 * Volume flux kernel for the six Cartesian stress components in the
 * three-dimensional first-order elastic (velocity-stress) system. The parameter
 * "component" selects which stress equation (sxx, syy, szz, sxy, sxz, syz) to
 * assemble. Each component contributes grad(test) dot F where F is built from
 * the coupled velocity gradients and Lame parameters.
 */
class ElasticStressFlux3D : public Kernel
{
public:
  static InputParameters validParams();
  ElasticStressFlux3D(const InputParameters & p);
protected:
  Real computeQpResidual() override;
  Real computeQpJacobian() override;
  Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _ux; const VariableValue & _uy; const VariableValue & _uz;
  const MaterialProperty<Real> & _lambda; const MaterialProperty<Real> & _mu;
  const unsigned int _comp; // 0..5
  // Coupled variable numbers for off-diagonal mapping
  const unsigned int _ux_var;
  const unsigned int _uy_var;
  const unsigned int _uz_var;
};
