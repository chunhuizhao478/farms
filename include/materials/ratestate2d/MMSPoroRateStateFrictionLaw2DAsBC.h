#pragma once

#include "PoroCZMComputeLocalTractionTotalBaseRSF2D.h"

class MMSPoroRateStateFrictionLaw2DAsBC : public PoroCZMComputeLocalTractionTotalBaseRSF2D
{
public:
  static InputParameters validParams();
  MMSPoroRateStateFrictionLaw2DAsBC(const InputParameters & parameters);

protected:
  /// method computing the total traction and its derivatives
  void computeInterfaceTractionAndDerivatives() override;

  //element length
  Real _len;

  //rate-and-state friction coefficients
  Real _f_o;
  Real _rsf_a;
  Real _rsf_b;
  Real _rsf_L;
  Real _delta_o;

  //density
  const MaterialProperty<Real> & _density;
  const MaterialProperty<Real> & _rhof;

  const VariableValue & _len2;

  //rotation matrix
  const MaterialProperty<RankTwoTensor> & _rot;
   
  //restoration forces
  const VariableValue & _reaction_rsf_x;
  const VariableValue & _reaction_rsf_y;
  const VariableValue & _reaction_rsf_neighbor_x;
  const VariableValue & _reaction_rsf_neighbor_y;

  //restoration pressures to calculate effective stress
  const VariableValue & _reaction_rsf_pressure_x;
  const VariableValue & _reaction_rsf_pressure_y;
  const VariableValue & _reaction_rsf_neighbor_pressure_x;
  const VariableValue & _reaction_rsf_neighbor_pressure_y;

  //restoration forces
  const VariableValue & _reaction_damp_x;
  const VariableValue & _reaction_damp_y;
  const VariableValue & _reaction_damp_neighbor_x;
  const VariableValue & _reaction_damp_neighbor_y;

  const VariableValue & _reaction_pressdamp_x;
  const VariableValue & _reaction_pressdamp_y;
  const VariableValue & _reaction_pressdamp_neighbor_x;
  const VariableValue & _reaction_pressdamp_neighbor_y;

    // Permeability and interface pressure
  const VariableValue & _interface_pressure_plus;
  const VariableValue & _interface_pressure_minus;

  

  //shear stress perturbatio

  //old along fault properties
  const MaterialProperty<Real> & _alongfaultvel_strike_plus_old;
  const MaterialProperty<Real> & _alongfaultvel_strike_minus_old;
  const MaterialProperty<Real> & _alongfaultvel_normal_plus_old;
  const MaterialProperty<Real> & _alongfaultvel_normal_minus_old;
  const MaterialProperty<Real> & _alongfaultdisp_strike_plus_old;
  const MaterialProperty<Real> & _alongfaultdisp_strike_minus_old;
  const MaterialProperty<Real> & _alongfaultdisp_normal_plus_old;
  const MaterialProperty<Real> & _alongfaultdisp_normal_minus_old;
  
  const MaterialProperty<Real> & _alongfaultdisp_strike_plus_older;
  const MaterialProperty<Real> & _alongfaultdisp_normal_plus_older;
  const MaterialProperty<Real> & _alongfaultdisp_strike_minus_older;
  const MaterialProperty<Real> & _alongfaultdisp_normal_minus_older;

  const MaterialProperty<Real> & _sliprate_strike_old;
  const MaterialProperty<Real> & _sliprate_normal_old;
  const MaterialProperty<Real> & _sliprate_mag_old;

  const MaterialProperty<Real> & _sliprate_predict;

  const MaterialProperty<Real> & _slip_strike_old;
  const MaterialProperty<Real> & _slip_normal_old;

  //old state variable
  const MaterialProperty<Real> & _statevar_old;
  const MaterialProperty<Real> & _statevar_older;
  
  //old traction
  const MaterialProperty<Real> & _traction_strike_old;
  const MaterialProperty<Real> & _traction_normal_old;


  const VariableValue & _mms_sliprate;
   const VariableValue & _mms_shear;
   const VariableValue & _mms_normal;
   const VariableValue & _mms_pressure;
   

  // Fluid displacement and velocity
  const VariableValue & _fluid_disp_x;
  const VariableValue & _fluid_disp_neighbor_x;
  const VariableValue & _fluid_disp_y;
  const VariableValue & _fluid_disp_neighbor_y;
  const VariableValue & _fluid_vel_x;
  const VariableValue & _fluid_vel_neighbor_x;
  const VariableValue & _fluid_vel_y;
  const VariableValue & _fluid_vel_neighbor_y;

  MaterialProperty<Real> & _traction_strike_TSN;
  MaterialProperty<Real> & _traction_normal_TSN;
  
};