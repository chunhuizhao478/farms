#pragma once

#include "Material.h"

class StabilizationMaterial : public Material
{
public:
  static InputParameters validParams();
  
  StabilizationMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

private:
  MaterialProperty<Real> & _tau_pspg;
  const MaterialProperty<Real> & _shear_modulus;  // From your damage material
  const Real _coeff;
};