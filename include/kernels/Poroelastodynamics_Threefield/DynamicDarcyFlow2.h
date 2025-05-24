#pragma once

#include "TimeKernel.h"
#include "Material.h"

//Forward Declarations
class TimeIntegrator;
//Forward Declarations
class DynamicDarcyFlow2 : public TimeKernel
{
public:
  static InputParameters validParams();

  DynamicDarcyFlow2(const InputParameters & parameters);
  
  virtual void computeJacobian() override;

protected:
virtual Real computeQpResidual() override;
virtual Real computeQpJacobian() override;
virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
private:
bool _lumping;
const MaterialProperty<Real> & _rhof; // fluid density
const MaterialProperty<Real> & _taut; // fluid density
const MaterialProperty<Real> & _nf; // porosity
const MaterialProperty<Real> & _K; // hydraulic conductivity
unsigned int _us_var_num; // id of skeleton displacement variable
const VariableValue * _us_dot_dot; // skeleton acceleration for explicit
const VariableValue * _dus_dot_dot_du;
const VariableValue * _u_dot_factor;
const VariableValue * _u_dotdot_factor;
const VariableValue * _u_dot_old;
const VariableValue * _du_dot_du;
const VariableValue * _du_dotdot_du;
const VariableValue * _u_dot_factor_dof;
/// The TimeIntegrator
const TimeIntegrator & _time_integrator;
};
