#pragma once

#include "Kernel.h"

class PoroStabilization : public Kernel
{
public:
  static InputParameters validParams();
  
  PoroStabilization(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  // Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;
  
  // Material properties from your custom material
  const MaterialProperty<Real> & _tau_pspg;
  const MaterialProperty<RankTwoTensor> & _stress;  // Total stress
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;  // YOUR tangent modulus
  const MaterialProperty<RankTwoTensor> & _stress_off_diag_jacobian;  // ∂σ/∂p from YOUR material
};