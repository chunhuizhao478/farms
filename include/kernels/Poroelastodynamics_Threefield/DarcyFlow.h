#pragma once

#include "TimeKernel.h"
#include "Material.h"

class TimeIntegrator;

//Forward Declarations
class DarcyFlow : public TimeKernel
{
public:
  static InputParameters validParams();

  DarcyFlow(const InputParameters & parameters);

protected:
virtual Real computeQpResidual() override;
virtual Real computeQpJacobian() override;
private:
const MaterialProperty<Real> & _K; // hydraulic conductivity
};
