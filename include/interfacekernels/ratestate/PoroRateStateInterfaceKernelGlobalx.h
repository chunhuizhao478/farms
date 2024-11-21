#pragma once
#include "InterfaceKernel.h"
#include "JvarMapInterface.h"
/**
 * Rate-and-State Friction Interface Displacement Transfer
 */
class PoroRateStateInterfaceKernelGlobalx : public InterfaceKernel
{
public:
  static InputParameters validParams();
  PoroRateStateInterfaceKernelGlobalx(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  const std::string _base_name;

  // Primary variables
  const MaterialProperty<Real> & _disp_strike_plus;
  const MaterialProperty<Real> & _disp_strike_minus;

  // Derivatives with respect to displacements
  const MaterialProperty<Real> & _dalongfaultdisp_strike_plus_tplusdt_du_strike_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_plus_tplusdt_du_strike_minus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_plus_tplusdt_du_normal_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_plus_tplusdt_du_normal_minus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_minus_tplusdt_du_strike_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_minus_tplusdt_du_strike_minus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_minus_tplusdt_du_normal_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_minus_tplusdt_du_normal_minus;

  // Derivatives with respect to fluid velocities
  const MaterialProperty<Real> & _dalongfaultdisp_strike_plus_tplusdt_dvf_strike_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_plus_tplusdt_dvf_strike_minus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_plus_tplusdt_dvf_normal_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_plus_tplusdt_dvf_normal_minus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_minus_tplusdt_dvf_strike_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_minus_tplusdt_dvf_strike_minus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_minus_tplusdt_dvf_normal_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_minus_tplusdt_dvf_normal_minus;

  // Derivatives with respect to pressure
  const MaterialProperty<Real> & _dalongfaultdisp_strike_plus_tplusdt_dp_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_plus_tplusdt_dp_minus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_minus_tplusdt_dp_plus;
  const MaterialProperty<Real> & _dalongfaultdisp_strike_minus_tplusdt_dp_minus;

  // Coupled variable numbers
  const unsigned int _y_var;
  const unsigned int _vfx_var;
  const unsigned int _vfy_var;
  const unsigned int _pressure_var;

  // Variable data for accessing shape functions
  MooseVariable & _vfx_var_data;
  MooseVariable & _vfy_var_data;
  MooseVariable & _pressure_var_data;

  // Shape functions for different variables
  const VariablePhiValue & _phi_vfx;
  const VariablePhiValue & _phi_vfx_neighbor;
  const VariablePhiValue & _phi_vfy;
  const VariablePhiValue & _phi_vfy_neighbor;
  const VariablePhiValue & _phi_pressure;
  const VariablePhiValue & _phi_pressure_neighbor;
};