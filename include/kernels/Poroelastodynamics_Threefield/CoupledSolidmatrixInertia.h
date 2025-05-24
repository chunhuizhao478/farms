#pragma once

#include "TimeKernel.h"
#include "Material.h"
#include "InertialForce.h"

//Forward Declarations
class TimeIntegrator;
//Forward Declarations
class CoupledSolidmatrixInertia : public TimeKernel
{
public:
  static InputParameters validParams();

  CoupledSolidmatrixInertia(const InputParameters & parameters);

protected:
virtual Real computeQpResidual() override;
virtual Real computeQpJacobian() override;
virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
private:
bool _lumping;
const MaterialProperty<Real> & _rhof; // fluid density
unsigned int _us_var_num; // id of skeleton displacement variable
const VariableValue & _us_older;
const VariableValue & _us_old;
const VariableValue & _us;
const VariableValue * _us_dot_dot; // skeleton acceleration for explicit
const VariableValue * _us_dot_dot_dof; // skeleton acceleration for explicit
const VariableValue * _dus_dot_dot_du;
const VariableValue * _u_dot_factor;
const VariableValue * _du_dot_du;
const TimeIntegrator & _time_integrator;;
};
