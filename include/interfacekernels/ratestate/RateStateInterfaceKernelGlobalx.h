#pragma once

#include "InterfaceKernel.h"

/**
 * Rate-and-State Friction Interface Displacement Transfer
 */
class RateStateInterfaceKernelGlobalx : public InterfaceKernel
{
public:
  static InputParameters validParams();

  RateStateInterfaceKernelGlobalx(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  const MaterialProperty<Real> & _disp_strike_plus;
  const MaterialProperty<Real> & _disp_strike_minus;
  
  // Add the derivatives
  const MaterialProperty<Real> & _dalongfaultdisp_strike_plus_tplusdt_du_strike_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_plus_tplusdt_du_strike_minus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_plus_tplusdt_du_normal_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_plus_tplusdt_du_normal_minus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_minus_tplusdt_du_strike_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_minus_tplusdt_du_strike_minus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_minus_tplusdt_du_normal_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_minus_tplusdt_du_normal_minus;

  // Coupled variable number for y-direction
  const unsigned int _y_var;
};