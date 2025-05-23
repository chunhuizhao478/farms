// CoupledFluidInertialForce.h
#pragma once
#include "TimeKernel.h"
#include "Material.h"

class TimeIntegrator;

class CoupledFluidInertialForce : public TimeKernel
{
public:
  static InputParameters validParams();
  CoupledFluidInertialForce(const InputParameters & parameters);

protected:
  bool _lumping;
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  const MaterialProperty<Real> & _rhof; // fluid density
  
  // Variable value pointers - correct declaration
  const VariableValue & _wf;           // Current fluid velocity
  const VariableValue & _wf_older;     // Old fluid velocity
  
  // Vector tag values for time integration
  const VariableValue * _wf_dot_factor;      // fluid relative acceleration
  const VariableValue * _wf_dot_factor_dof;  // fluid relative acceleration DOF
  const VariableValue * _dwf_dot_du;         // fluid relative acceleration derivative
  
  unsigned int _w_var_num;              // id of the Darcy vel variable
  
  /// The TimeIntegrator
  TimeIntegrator & _time_integrator;
};