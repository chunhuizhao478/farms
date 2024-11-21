#pragma once

#include "InterfaceKernel.h"

/**
 * Rate-and-State Friction Interface Displacement Transfer
 */
class RateStateInterfaceKernelGlobaly : public InterfaceKernel
{
public:
  static InputParameters validParams();

  RateStateInterfaceKernelGlobaly(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  const MaterialProperty<Real> & _disp_normal_plus;
  const MaterialProperty<Real> & _disp_normal_minus;
  
  // Add the derivatives
  const MaterialProperty<Real> & _dalongfaultdisp_normal_plus_tplusdt_du_strike_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_normal_plus_tplusdt_du_strike_minus;
  const MaterialProperty<Real> & _dalongfaultdisp_normal_plus_tplusdt_du_normal_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_normal_plus_tplusdt_du_normal_minus;
  const MaterialProperty<Real> & _dalongfaultdisp_normal_minus_tplusdt_du_strike_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_normal_minus_tplusdt_du_strike_minus;
  const MaterialProperty<Real> & _dalongfaultdisp_normal_minus_tplusdt_du_normal_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_normal_minus_tplusdt_du_normal_minus;
  
  // Coupled variable number for x-direction
  const unsigned int _x_var;
};