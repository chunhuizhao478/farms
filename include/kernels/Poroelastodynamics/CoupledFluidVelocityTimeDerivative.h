#pragma once

#include "TimeKernel.h"
#include "PorousFlowDictator.h"
#include "RankTwoTensor.h"

/**
 * This calculates the time derivative for a coupled variable
 **/
class CoupledFluidVelocityTimeDerivative : public TimeKernel
{
public:
  static InputParameters validParams();

  CoupledFluidVelocityTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  const VariableValue & _v_dot;
  const VariableValue & _dv_dot;
  const unsigned int _v_var;

  /// The fluid component index
  const unsigned int _fluid_component;

  /// PorousFlowDictator UserObject
  const PorousFlowDictator & _dictator;

  /// Number of fluid phases
  const unsigned int _num_phases;

  /// Fluid density for each phase (at the qp)
  const MaterialProperty<std::vector<Real>> & _fluid_density_qp;

};