#pragma once

#include "AuxKernel.h"
//Forward Declarations
/**
* Accumulate values from one auxiliary variable into another
*/

class NewmarkPoreFluidAccelAux : public AuxKernel
{
public:
  static InputParameters validParams();
  NewmarkPoreFluidAccelAux(const InputParameters & parameters);
  
protected:
virtual Real computeValue();
const VariableValue & _w_old; // W is Darcy velocity
const VariableValue & _w;
const VariableValue & _u_old;
Real _gamma;
};
