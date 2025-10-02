#pragma once

#include "Material.h"

class Function;

/**
 * FarmsFunctionMaterial evaluates a MOOSE Function at each quadrature point
 * and stores the result in a scalar material property.
 */
class FarmsFunctionMaterial : public Material
{
public:
  static InputParameters validParams();

  FarmsFunctionMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  const Function & _function;
  MaterialProperty<Real> & _property;
};

