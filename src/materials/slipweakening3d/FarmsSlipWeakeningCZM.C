//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FarmsSlipWeakeningCZM.h"
registerMooseObject("farmsApp", FarmsSlipWeakeningCZM);

InputParameters
FarmsSlipWeakeningCZM::validParams()
{
  InputParameters params = FarmsSlipWeakeningBase::validParams();
  params.addClassDescription("Linear mixed extrinsic cohesive law of "
                             "Ortiz-Pandolfi model");
  params.addParam<Real>("Dc", 1.0, "Value of characteristic length");
  params.addRequiredCoupledVar("disp_slipweakening_x","displacement in x dir");
  params.addRequiredCoupledVar("disp_slipweakening_y","displacement in y dir");
  params.addRequiredCoupledVar("disp_slipweakening_z","displacement in z dir");
  params.addRequiredCoupledVar("reaction_slipweakening_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_slipweakening_y","reaction in y dir");
  params.addRequiredCoupledVar("reaction_slipweakening_z","reaction in z dir");
  params.addRequiredCoupledVar("reaction_damp_x","reaction damping in x dir");
  params.addRequiredCoupledVar("reaction_damp_y","reaction damping in y dir");
  params.addRequiredCoupledVar("reaction_damp_z","reaction damping in z dir");
  params.addRequiredCoupledVar("elem_length","length of a element");
  params.addRequiredCoupledVar("mu_s","static friction coefficient spatial distribution");
  params.addRequiredCoupledVar("mu_d","dynamic friction coefficient spatial distribution");
  params.addCoupledVar("cohesion","cohesion in shear stress");
  params.addCoupledVar("forced_rupture_time","time of forced rupture");  
  return params;
}

FarmsSlipWeakeningCZM::FarmsSlipWeakeningCZM(const InputParameters & parameters)
  : FarmsSlipWeakeningBase(parameters),
  _Dc(getParam<Real>("Dc")),
  _density(getMaterialPropertyByName<Real>("density")),
  _rot(getMaterialPropertyByName<RealTensorValue>("rotation_matrix")),
  _disp_slipweakening_x(coupledValue("disp_slipweakening_x")),
  _disp_slipweakening_neighbor_x(coupledNeighborValue("disp_slipweakening_x")),
  _disp_slipweakening_y(coupledValue("disp_slipweakening_y")),
  _disp_slipweakening_neighbor_y(coupledNeighborValue("disp_slipweakening_y")),
  _disp_slipweakening_z(coupledValue("disp_slipweakening_z")),
  _disp_slipweakening_neighbor_z(coupledNeighborValue("disp_slipweakening_z")),
  _reaction_slipweakening_x(coupledValue("reaction_slipweakening_x")),
  _reaction_slipweakening_neighbor_x(coupledNeighborValue("reaction_slipweakening_x")),
  _reaction_slipweakening_y(coupledValue("reaction_slipweakening_y")),
  _reaction_slipweakening_neighbor_y(coupledNeighborValue("reaction_slipweakening_y")),
  _reaction_slipweakening_z(coupledValue("reaction_slipweakening_z")),
  _reaction_slipweakening_neighbor_z(coupledNeighborValue("reaction_slipweakening_z")),
  _reaction_damp_x(coupledValue("reaction_damp_x")),
  _reaction_damp_neighbor_x(coupledNeighborValue("reaction_damp_x")),
  _reaction_damp_y(coupledValue("reaction_damp_y")),
  _reaction_damp_neighbor_y(coupledNeighborValue("reaction_damp_y")),
  _reaction_damp_z(coupledValue("reaction_damp_z")),
  _reaction_damp_neighbor_z(coupledNeighborValue("reaction_damp_z")),
  _disp_slipweakening_x_old(coupledValueOld("disp_slipweakening_x")),
  _disp_slipweakening_neighbor_x_old(coupledNeighborValueOld("disp_slipweakening_x")),
  _disp_slipweakening_y_old(coupledValueOld("disp_slipweakening_y")),
  _disp_slipweakening_neighbor_y_old(coupledNeighborValueOld("disp_slipweakening_y")),
  _disp_slipweakening_z_old(coupledValueOld("disp_slipweakening_z")),
  _disp_slipweakening_neighbor_z_old(coupledNeighborValueOld("disp_slipweakening_z")),
  _elem_length(coupledValue("elem_length")),
  _mu_s(coupledValue("mu_s")),
  _mu_d(coupledValue("mu_d")),
  _sts_init(getMaterialPropertyByName<RankTwoTensor>("static_initial_stress_tensor_slipweakening")),
  _Co_coupled(isCoupled("cohesion")),
  _T_coupled(isCoupled("forced_rupture_time")),
  _Co(_Co_coupled ? &coupledValue("cohesion") : nullptr),
  _T(_T_coupled ? &coupledValue("forced_rupture_time") : nullptr),
  _accumulated_slip_along_normal_old(getMaterialPropertyOldByName<Real>("accumulated_slip_along_normal")),
  _accumulated_slip_along_strike_old(getMaterialPropertyOldByName<Real>("accumulated_slip_along_strike")),
  _accumulated_slip_along_dip_old(getMaterialPropertyOldByName<Real>("accumulated_slip_along_dip")),
  _slip_along_normal_old(getMaterialPropertyOldByName<Real>("slip_along_normal")),
  _slip_along_strike_old(getMaterialPropertyOldByName<Real>("slip_along_strike")),
  _slip_along_dip_old(getMaterialPropertyOldByName<Real>("slip_along_dip"))
{
}

void
FarmsSlipWeakeningCZM::initQpStatefulProperties()
{}

RealVectorValue
FarmsSlipWeakeningCZM::computeTraction()
{ 
  //Global Displacement Jump
  RealVectorValue displacement_jump_global(_disp_slipweakening_x[_qp]-_disp_slipweakening_neighbor_x[_qp],_disp_slipweakening_y[_qp]-_disp_slipweakening_neighbor_y[_qp],_disp_slipweakening_z[_qp]-_disp_slipweakening_neighbor_z[_qp]);
  _displacement_jump_global[_qp] = displacement_jump_global;

  //Global Displacement Jump Old
  RealVectorValue displacement_jump_old_global(_disp_slipweakening_x_old[_qp]-_disp_slipweakening_neighbor_x_old[_qp],_disp_slipweakening_y_old[_qp]-_disp_slipweakening_neighbor_y_old[_qp],_disp_slipweakening_z_old[_qp]-_disp_slipweakening_neighbor_z_old[_qp]);

  //Global Displacement Jump Rate
  RealVectorValue displacement_jump_rate_global = (displacement_jump_global - displacement_jump_old_global)*(1/_dt);  
  _displacement_jump_rate_global[_qp] = displacement_jump_rate_global;

  //Local Displacement Jump / Displacement Jump Rate
  RealVectorValue displacement_jump_local = GlobaltoLocalVector(displacement_jump_global, _rot[_qp]);
  RealVectorValue displacement_jump_rate_local = GlobaltoLocalVector(displacement_jump_rate_global, _rot[_qp]);

  //Parameter Initialization
  Real mu_s = _mu_s[_qp]; 
  Real mu_d = _mu_d[_qp]; 
  Real Dc = _Dc; 
  Real tau_f = 0;  

  //Background Stress Projection
  RealTensorValue sts_init_global = RankTwoTensor2RealTensorValue(_sts_init[_qp]);
  RealTensorValue sts_init_local = GlobaltoLocalMatrix(sts_init_global, _rot[_qp]);

  //Local Initial Traction
  RealVectorValue traction_init_local = MatrixVectorMultiply(sts_init_local, _normals[_qp]);

  Real T_strike_o = traction_init_local(0);
  Real T_dip_o    = traction_init_local(1);
  Real T_normal_o = traction_init_local(2);

  //*Restoration Force*
  //Stress Divergence Components (label as stsdivcomp)
  //--------------------------------------------------------------------------------------------------//  

  ///Define in global coordinate
  //current time step 
  RealVectorValue R_plus_global_stsdivcomp(-_reaction_slipweakening_x[_qp],-_reaction_slipweakening_y[_qp], -_reaction_slipweakening_z[_qp]);
  RealVectorValue R_minus_global_stsdivcomp(-_reaction_slipweakening_neighbor_x[_qp],-_reaction_slipweakening_neighbor_y[_qp], -_reaction_slipweakening_neighbor_z[_qp]);

  ///Rotate in local coordinate
  //current time step
  RealVectorValue R_plus_local_stsdivcomp  = GlobaltoLocalVector(R_plus_global_stsdivcomp, _rot[_qp]);
  RealVectorValue R_minus_local_stsdivcomp = GlobaltoLocalVector(R_minus_global_stsdivcomp, _rot[_qp]);

  ///Get Components
  //current time step  
  Real R_plus_local_strike_stsdivcomp  = R_plus_local_stsdivcomp(0);
  Real R_plus_local_dip_stsdivcomp     = R_plus_local_stsdivcomp(1);
  Real R_plus_local_normal_stsdivcomp  = R_plus_local_stsdivcomp(2);

  Real R_minus_local_strike_stsdivcomp = R_minus_local_stsdivcomp(0);
  Real R_minus_local_dip_stsdivcomp    = R_minus_local_stsdivcomp(1);
  Real R_minus_local_normal_stsdivcomp = R_minus_local_stsdivcomp(2);

  //--------------------------------------------------------------------------------------------------//

  //Damping Components Contribution (label as dampingcomp)

  ///Define in global coordinate
  //current time step 
  RealVectorValue R_plus_global_dampingcomp(-_reaction_damp_x[_qp],-_reaction_damp_y[_qp], -_reaction_damp_z[_qp]);
  RealVectorValue R_minus_global_dampingcomp(-_reaction_damp_neighbor_x[_qp],-_reaction_damp_neighbor_y[_qp], -_reaction_damp_neighbor_z[_qp]);  

  ///Rotate in local coordinate
  //current time step
  RealVectorValue R_plus_local_dampingcomp = GlobaltoLocalVector(R_plus_global_dampingcomp, _rot[_qp]);
  RealVectorValue R_minus_local_dampingcomp = GlobaltoLocalVector(R_minus_global_dampingcomp, _rot[_qp]);

  ///Get Components
  //current time step  
  Real R_plus_local_strike_dampingcomp  = R_plus_local_dampingcomp(0);
  Real R_plus_local_dip_dampingcomp     = R_plus_local_dampingcomp(1);
  Real R_plus_local_normal_dampingcomp  = R_plus_local_dampingcomp(2);

  Real R_minus_local_strike_dampingcomp = R_minus_local_dampingcomp(0);
  Real R_minus_local_dip_dampingcomp    = R_minus_local_dampingcomp(1);
  Real R_minus_local_normal_dampingcomp = R_minus_local_dampingcomp(2);

  //--------------------------------------------------------------------------------------------------//

  //Reaction force in local coordinate
  Real R_plus_local_strike  = R_plus_local_strike_stsdivcomp  + R_plus_local_strike_dampingcomp;
  Real R_plus_local_dip     = R_plus_local_dip_stsdivcomp     + R_plus_local_dip_dampingcomp;
  Real R_plus_local_normal  = R_plus_local_normal_stsdivcomp  + R_plus_local_normal_dampingcomp;
  Real R_minus_local_strike = R_minus_local_strike_stsdivcomp + R_minus_local_strike_dampingcomp;
  Real R_minus_local_dip    = R_minus_local_dip_stsdivcomp    + R_minus_local_dip_dampingcomp;  
  Real R_minus_local_normal = R_minus_local_normal_stsdivcomp + R_minus_local_normal_dampingcomp;
  
  //element length
  Real elem_length = _elem_length[_qp];

  //Compute node mass //equal length tetrahedron
  Real M = _density[_qp] * elem_length * elem_length * elem_length / 8;

  //Compute sticking stress
  Real Tstrike = (1/_dt)*M*displacement_jump_rate_local(0)/(2*elem_length*elem_length) + (R_plus_local_strike - R_minus_local_strike)/(2*elem_length*elem_length) + T_strike_o;
  Real Tdip    = (1/_dt)*M*displacement_jump_rate_local(1)/(2*elem_length*elem_length) + (R_plus_local_dip    - R_minus_local_dip   )/(2*elem_length*elem_length) + T_dip_o;
  Real Tnormal = (1/_dt)*M*(displacement_jump_rate_local(2)+(1/_dt)*displacement_jump_local(2))/(2*elem_length*elem_length) + ( (R_plus_local_normal - R_minus_local_normal) / ( 2*elem_length*elem_length ) ) + T_normal_o ;  

  //Compute fault traction
  //min(0,sigma_N)
  if (Tnormal<0)
  {
  }else{
    Tnormal = 0;
  }

  //Compute friction strength
  if (std::norm(displacement_jump_local) < Dc)
  {
    tau_f = (mu_s - (mu_s - mu_d)*std::norm(displacement_jump_local)/Dc)*(-Tnormal); // square for shear component
  }
  else
  {
    tau_f = mu_d * (-Tnormal);
  }  

  //Compute fault traction
  if (std::sqrt(Tstrike*Tstrike + Tdip*Tdip)<tau_f)
  {
  }else{
    Tstrike = tau_f*Tstrike/std::sqrt(Tstrike*Tstrike + Tdip*Tdip);
    Tdip    = tau_f*Tdip/std::sqrt(Tstrike*Tstrike + Tdip*Tdip);
  }  

  //Assign back traction in CZM
  RealVectorValue traction_local(0.0);

  traction_local(0) = Tstrike - T_strike_o; 
  traction_local(1) = Tdip    - T_dip_o; 
  traction_local(2) = Tnormal - T_normal_o;

  //Rotate back traction difference to global coordinates
  RealVectorValue traction_global(0.0);
  traction_global = LocaltoGlobalVector(traction_local, _rot[_qp]);

  return traction_global;
}

RealTensorValue
FarmsSlipWeakeningCZM::computeTractionDerivatives()
{ 
  RealTensorValue tangent_modulus_on_interface(0.0);
  return tangent_modulus_on_interface;
}
