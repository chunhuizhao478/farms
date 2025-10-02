#pragma once
#include "Kernel.h"

/**
 * Volume flux kernel for the velocity equations of the three-dimensional
 * first-order elastic system. The parameter "component" chooses ux, uy, or uz
 * and the kernel assembles grad(test) dot the matching stress column.
 */
class ElasticVelocityFlux3D : public Kernel
{
public:
  static InputParameters validParams();
  ElasticVelocityFlux3D(const InputParameters & p);
protected:
  Real computeQpResidual() override;
  // Diagonal Jacobian is zero (flux does not depend on the velocity unknown itself)
  Real computeQpJacobian() override { return 0.0; }
  // Off-diagonal derivatives wrt stress components
  Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // Coupled stress tensor components
  const VariableValue & _sxx; const VariableValue & _syy; const VariableValue & _szz;
  const VariableValue & _sxy; const VariableValue & _sxz; const VariableValue & _syz;
  const unsigned int _comp; // 0=u,1=v,2=w
  // Variable numbers (initialized in constructor)
  const unsigned int _sxx_var;
  const unsigned int _syy_var;
  const unsigned int _szz_var;
  const unsigned int _sxy_var;
  const unsigned int _sxz_var;
  const unsigned int _syz_var;
};
