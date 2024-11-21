#pragma once

#include "Kernel.h"
#include "Material.h"
//Forward Declarations
class DynamicDarcyFlow : public Kernel
{
public:
  static InputParameters validParams();

  DynamicDarcyFlow(const InputParameters & parameters);
  
protected:
virtual Real computeQpResidual();
virtual Real computeQpJacobian();
virtual Real computeQpOffDiagJacobian(unsigned int jvar);
private:
const MaterialProperty<Real> & _rhof; // fluid density
const MaterialProperty<Real> & _nf; // porosity
const MaterialProperty<Real> & _K; // hydraulic conductivity
const VariableValue & _us; // skeleton displacement
const VariableValue & _us_old;
const VariableValue & _vs_old; // skeleton velocity
const VariableValue & _as_old; // skeleton acceleration
const VariableValue & _u_old; // Darcy velocity
// this is actually w, but is called u in this kernel
// because this the variable for this kernel
const VariableValue & _af_old; // fluid relative acceleration
unsigned int _us_var_num; // id of skeleton displacement variable
const Real _gravity; // acceleration due to gravity
const Real _beta;
const Real _gamma;
};
