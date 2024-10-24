#pragma once

#include "CZMComputeLocalTractionTotalBaseRSF2D.h"

class RateStateFrictionLaw2Dv7 : public CZMComputeLocalTractionTotalBaseRSF2D
{
public:
  static InputParameters validParams();
  RateStateFrictionLaw2Dv7(const InputParameters & parameters);

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

  //rotation matrix
  const MaterialProperty<RankTwoTensor> & _rot;
   
  //restoration forces
  const VariableValue & _reaction_rsf_x;
  const VariableValue & _reaction_rsf_y;
  const VariableValue & _reaction_rsf_neighbor_x;
  const VariableValue & _reaction_rsf_neighbor_y;

  //shear stress perturbation
  ///Measure from current time step
  const VariableValue & _Ts_perturb;
  const VariableValue & _Ts_perturb_old;

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
  const MaterialProperty<Real> & _alongfaultdisp_strike_minus_older;
  const MaterialProperty<Real> & _alongfaultdisp_normal_plus_older;
  const MaterialProperty<Real> & _alongfaultdisp_normal_minus_older;
  
  const MaterialProperty<Real> & _sliprate_strike_old;
  const MaterialProperty<Real> & _sliprate_normal_old;
  const MaterialProperty<Real> & _sliprate_mag_old;
  
  const MaterialProperty<Real> & _slip_strike_old;
  const MaterialProperty<Real> & _slip_normal_old;
  const MaterialProperty<Real> & _slip_mag_old;
  
  //old state variable
  const MaterialProperty<Real> & _statevar_old;
  
  //old traction
  const MaterialProperty<Real> & _traction_strike_old;
  const MaterialProperty<Real> & _traction_normal_old;
};