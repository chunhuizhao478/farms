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
const VariableValue & _wf_older;
const VariableValue & _wf;
const VariableValue * _wf_dot_factor; // fluid relative acceleration
const VariableValue * _wf_dot_factor_dof; // fluid relative acceleration
const VariableValue * _dwf_dot_du; // fluid relative acceleration
unsigned int _w_var_num; // id of the Darcy vel variable
/// The TimeIntegrator
TimeIntegrator & _time_integrator;
};