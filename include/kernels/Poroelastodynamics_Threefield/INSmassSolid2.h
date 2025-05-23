#pragma once

#include "Kernel.h"
//Forward Declarations
class INSmassSolid2 : public Kernel
{
public:
  static InputParameters validParams();

  INSmassSolid2(const InputParameters & parameters);
protected:
virtual Real computeQpResidual();
virtual Real computeQpJacobian();
virtual Real computeQpOffDiagJacobian(unsigned int jvar);
private:

const VariableValue & _u_dot; 
const VariableValue & _du_dot_du; 
const MaterialProperty<Real> & _coefficient_M; 
unsigned int _ndisp; // number of displacement components (1D, 2D or 3D)
unsigned int _ux_var; // id of displacement components
unsigned int _uy_var;
unsigned int _uz_var;
const VariableGradient &_grad_ux; // gradient of displacement
const VariableGradient &_grad_uy;
const VariableGradient &_grad_uz;
const VariableGradient &_grad_ux_older; // gradient of previous displacement
const VariableGradient &_grad_uy_older;
const VariableGradient &_grad_uz_older;
const MaterialProperty<Real> & _coefficient;  /// Biot coefficient
};