#pragma once

#include "Kernel.h"
#include "Material.h"
//Forward Declarations
class FluidStorageNewmark : public Kernel
{
public:
  static InputParameters validParams();

  FluidStorageNewmark(const InputParameters & parameters);
  
protected:
virtual Real computeQpResidual();
virtual Real computeQpJacobian();
virtual Real computeQpOffDiagJacobian(unsigned int jvar);
private:
const MaterialProperty<Real> & _coefficient; 
const VariableValue & _u_old; 
};
