#pragma once
#include "TimeKernel.h"
#include "Material.h"

class DynamicDarcyFlow3 : public TimeKernel
{
public:
  static InputParameters validParams();
  DynamicDarcyFlow3(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  // Material properties
  const MaterialProperty<Real> & _rhof;  // fluid density
  const MaterialProperty<Real> & _taut;  // tortuosity
  const MaterialProperty<Real> & _nf;    // porosity
  const MaterialProperty<Real> & _K;     // hydraulic conductivity

  // Current and previous time step values for fluid variable
  const VariableValue & _u_older;        // u at two timesteps ago

  // Current and previous time step values for skeleton acceleration
  const VariableValue & _us;             // skeleton acceleration at current timestep
  const VariableValue & _us_old;         // skeleton acceleration at previous timestep
  const VariableValue & _us_older;       // skeleton acceleration at two timesteps ago
  
  // Variable number for skeleton acceleration (needed for off-diagonal Jacobian)
  unsigned int _us_var_num;
};