#pragma once

#include "Kernel.h"
#include "Material.h"
//Forward Declarations
class PoreFluidInertialForceCoupling : public Kernel
{
public:
  static InputParameters validParams();

  PoreFluidInertialForceCoupling(const InputParameters & parameters);

protected:
virtual Real computeQpResidual();
virtual Real computeQpJacobian();
virtual Real computeQpOffDiagJacobian(unsigned int jvar);
private:
const MaterialProperty<Real> & _rhof; // fluid density
const VariableValue & _af_old; // previous value of fluid acceleration
const VariableValue & _w; // Darcy velocity
const VariableValue & _w_old;
unsigned int _w_var_num; // id of the Darcy vel variable
const Real _gamma;
};