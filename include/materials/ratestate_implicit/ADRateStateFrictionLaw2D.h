#pragma once

#include "ADCZMComputeLocalTractionTotalBaseRSF2D.h"

class ADRateStateFrictionLaw2D : public ADCZMComputeLocalTractionTotalBaseRSF2D
{
public:
  static InputParameters validParams();
  ADRateStateFrictionLaw2D(const InputParameters & parameters);

protected:
  virtual void computeInterfaceTraction();

  //rate-and-state friction coefficients
  ADReal _f_o;
  ADReal _rsf_a;
  ADReal _rsf_b;
  ADReal _rsf_L;
  ADReal _delta_o;
  ADReal _Vini;

  //initial shear and normal traction
  ADReal _Tn_o;
  ADReal _Ts_o;

  //old/older interface jump
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_old;
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_older;

  //old statevar
  const MaterialProperty<Real> & _statevar_old;
  
};