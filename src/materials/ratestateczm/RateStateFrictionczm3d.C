//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*
Material Description of Slip Weakening Friction 3d
*/

#include "RateStateFrictionczm3d.h"
#include "InterfaceKernel.h"

registerMooseObject("farmsApp", RateStateFrictionczm3d);

InputParameters
RateStateFrictionczm3d::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("linear slip weakening traction separation law.");
  params.addRequiredParam<Real>("T1_o", "background shear traction in strike dir");
  params.addRequiredParam<Real>("T2_o", "background normal traction");
  params.addRequiredParam<Real>("T3_o", "background shear traction in dip dir");
  params.addRequiredParam<Real>("len","element length");
  params.addRequiredParam<Real>("f_o","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("rsf_a","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("rsf_b","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("rsf_L","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("delta_o","slip rate parameter");
  params.addRequiredParam<Real>("statevar_init","initial value of state variable");
  params.addRequiredParam<Real>("sliprate_strike_init","initial value of strike slip rate");
  params.addRequiredCoupledVar("disp_x", "displacement in x dir");
  params.addRequiredCoupledVar("disp_y", "displacement in y dir");
  params.addRequiredCoupledVar("disp_z", "displacement in z dir");
  params.addRequiredCoupledVar("vel_x","velocity in x dir");
  params.addRequiredCoupledVar("vel_y","velocity in y dir");
  params.addRequiredCoupledVar("vel_z","velocity in z dir");
  params.addRequiredCoupledVar("reaction_x", "reaction in x dir");
  params.addRequiredCoupledVar("reaction_y", "reaction in y dir");
  params.addRequiredCoupledVar("reaction_z", "reaction in z dir");
  params.addRequiredCoupledVar("Ts_perturb","shear stress perturbation in strike dir at time t");
  return params;
}

RateStateFrictionczm3d::RateStateFrictionczm3d(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _statevar(declareProperty<Real>("statevar")),
    _statevar_old(getMaterialPropertyOldByName<Real>("statevar")),
    _slipratevar(declareProperty<RealVectorValue>("slipratevar")),
    _slipratevar_old(getMaterialPropertyOldByName<RealVectorValue>("slipratevar")),
    _T1_o(getParam<Real>("T1_o")),
    _T2_o(getParam<Real>("T2_o")),
    _T3_o(getParam<Real>("T3_o")),
    _len(getParam<Real>("len")),
    _f_o(getParam<Real>("f_o")),
    _rsf_a(getParam<Real>("rsf_a")),
    _rsf_b(getParam<Real>("rsf_b")),
    _rsf_L(getParam<Real>("rsf_L")),
    _delta_o(getParam<Real>("delta_o")),
    _statevar_init(getParam<Real>("statevar_init")),
    _sliprate_strike_init(getParam<Real>("sliprate_strike_init")),
    _density(getMaterialPropertyByName<Real>(_base_name + "density")),
    _rot(getMaterialPropertyByName<RankTwoTensor>(_base_name + "czm_total_rotation")),
    _disp_x(coupledValue("disp_x")),
    _disp_neighbor_x(coupledNeighborValue("disp_x")),
    _disp_y(coupledValue("disp_y")),
    _disp_neighbor_y(coupledNeighborValue("disp_y")),
    _disp_z(coupledValue("disp_z")),
    _disp_neighbor_z(coupledNeighborValue("disp_z")),
    _vel_x(coupledValue("vel_x")),
    _vel_neighbor_x(coupledNeighborValue("vel_x")),
    _vel_y(coupledValue("vel_y")),
    _vel_neighbor_y(coupledNeighborValue("vel_y")),
    _vel_z(coupledValue("vel_z")),
    _vel_neighbor_z(coupledNeighborValue("vel_z")),
    _reaction_x(coupledValue("reaction_x")),
    _reaction_neighbor_x(coupledNeighborValue("reaction_x")),
    _reaction_y(coupledValue("reaction_y")),
    _reaction_neighbor_y(coupledNeighborValue("reaction_y")),
    _reaction_z(coupledValue("reaction_z")),
    _reaction_neighbor_z(coupledNeighborValue("reaction_z")),
    _disp_x_old(coupledValueOld("disp_x")),
    _disp_neighbor_x_old(coupledNeighborValueOld("disp_x")),
    _disp_y_old(coupledValueOld("disp_y")),
    _disp_neighbor_y_old(coupledNeighborValueOld("disp_y")),
    _disp_z_old(coupledValueOld("disp_z")),
    _disp_neighbor_z_old(coupledNeighborValueOld("disp_z")),
    _Ts_perturb(coupledValue("Ts_perturb")),
    _Ts_perturb_old(coupledValueOld("Ts_perturb"))
{

  // only works for small strain
  if (hasBlockMaterialProperty<RankTwoTensor>(_base_name + "strain_increment"))
  {
    mooseError("RateState only works for small strain!");
  }
}

void
RateStateFrictionczm3d::initQpStatefulProperties()
{
  _statevar[_qp]           = _statevar_init;

  _slipratevar[_qp](0)     = 0;                     //normal
  _slipratevar[_qp](1)     = _sliprate_strike_init; //strike
  _slipratevar[_qp](2)     = 0;                     //dip

}

void
RateStateFrictionczm3d::computeInterfaceTractionAndDerivatives()
{
  // Global Displacement Jump
  RealVectorValue displacement_jump_global(
      _disp_x[_qp] - _disp_neighbor_x[_qp],
      _disp_y[_qp] - _disp_neighbor_y[_qp],
      _disp_z[_qp] - _disp_neighbor_z[_qp]);

  // Global Displacement Jump Old
  RealVectorValue displacement_jump_old_global(
    _disp_x_old[_qp] - _disp_neighbor_x_old[_qp],
    _disp_y_old[_qp] - _disp_neighbor_y_old[_qp],
    _disp_z_old[_qp] - _disp_neighbor_z_old[_qp]);

  // Global Displacement Jump Rate
  RealVectorValue displacement_jump_rate_global =
      (displacement_jump_global - displacement_jump_old_global) * (1 / _dt);

  // Local Displacement Jump / Displacement Jump Rate
  RealVectorValue displacement_jump = _rot[_qp].transpose() * displacement_jump_global;
  RealVectorValue displacement_jump_rate = _rot[_qp].transpose() * displacement_jump_rate_global;

  // n is along normal direction; t is along tangential direction; d is along dip direction
  Real displacement_jump_n = displacement_jump(0);
  //Real displacement_jump_t = displacement_jump(1);
  //Real displacement_jump_d = displacement_jump(2);
  Real displacement_jump_rate_n = displacement_jump_rate(0);
  Real displacement_jump_rate_t = displacement_jump_rate(1);
  Real displacement_jump_rate_d = displacement_jump_rate(2);

  // Reaction force in local coordinate
  RealVectorValue R_plus_global(-_reaction_x[_qp],
                                -_reaction_y[_qp],
                                -_reaction_z[_qp]);
  RealVectorValue R_minus_global(-_reaction_neighbor_x[_qp],
                                 -_reaction_neighbor_y[_qp],
                                 -_reaction_neighbor_z[_qp]);

  RealVectorValue R_plus_local = _rot[_qp].transpose() * R_plus_global;
  RealVectorValue R_minus_local = _rot[_qp].transpose() * R_minus_global;

  // n is along normal direction; t is along tangential direction; d is along dip direction
  Real R_plus_local_n = R_plus_local(0);
  Real R_plus_local_t = R_plus_local(1);
  Real R_plus_local_d = R_plus_local(2);
  Real R_minus_local_n = R_minus_local(0);
  Real R_minus_local_t = R_minus_local(1);
  Real R_minus_local_d = R_minus_local(2);

  // Compute node mass and area
  Real M = 0;
  Real A = 0;
  if (_current_elem->type() == libMesh::ElemType::TET4){
    M = (_density[_qp] * sqrt(2) * _len * _len * _len / 12 / 4) * 6;
    A = (sqrt(3) * _len * _len / 4 / 3) * 6;
  }
  else if (_current_elem->type() == libMesh::ElemType::HEX8){
    M = (_density[_qp] * _len * _len * _len / 8) * 4;
    A = (_len * _len / 4) * 4;
  }

  // Compute normal sticking stress
  Real T2 = -(1 / _dt) * M * (displacement_jump_rate_n + (1 / _dt) * displacement_jump_n) /
                (2 * A) +
            ((R_minus_local_n - R_plus_local_n) / (2 * A)) - _T2_o;

  // Compute fault traction
  if (T2 < 0)
  {
  }
  else
  {
    T2 = 0;
  }

  ///Make Tn positive
  Real Tn = abs(T2);

  //*Compute Trial Shear Traction Along Strike Direction at Current Time Step*
  Real Ts_trial = (1 / _dt) * M * displacement_jump_rate_t / (2 * A) +
                  (R_plus_local_t - R_minus_local_t) / (2 * A) + _T1_o + _Ts_perturb[_qp];
  Real Td_trial = (1 / _dt) * M * displacement_jump_rate_d / (2 * A) +
                  (R_plus_local_d - R_minus_local_d) / (2 * A) + _T3_o;
  Real Tmag_trial = sqrt(Ts_trial*Ts_trial+Td_trial*Td_trial);

  //*Compute slip rate magnitude old
  // _slipratevar (n,t,d)
  Real sliprate_mag_old = sqrt(_slipratevar_old[_qp](1)*_slipratevar_old[_qp](1)+_slipratevar_old[_qp](2)*_slipratevar_old[_qp](2));

  //const
  Real c = A * _dt * ( M + M ) / (M * M);
  Real Z = 0.5 / _delta_o * exp((_f_o + _rsf_b * log(_delta_o * _statevar_old[_qp]/_rsf_L))/_rsf_a);

  //Setup while loop
  Real iterr = 1;
  Real max_iter = 10000;
  Real er = 1;
  Real solution;
  Real guess_i = sliprate_mag_old; //slip rate at time t-dt/2
  Real residual;
  Real jacobian;
  Real guess_j;
  while ( er > 1e-10 && iterr < max_iter ){

      //Compute Residual
      residual = guess_i + c * Tn * _rsf_a * asinh( 0.5*(guess_i+sliprate_mag_old) * Z ) - c * Tmag_trial;

      //Compute Jacobian
      jacobian = 1.0 + c * Tn * _rsf_a * 0.5 * Z / sqrt( 1.0 + 0.5 * 0.5 * (guess_i+sliprate_mag_old) * (guess_i+sliprate_mag_old) * Z * Z );

      //Compute New guess
      guess_j = guess_i - residual / jacobian;

      //save
      solution = guess_j;

      //Compute err
      er = abs(guess_j - guess_i)/abs(guess_j);

      //Update Old guess
      guess_i = guess_j;

      //update iterr
      iterr = iterr + 1;

  }

  if (iterr == max_iter){
      mooseError("NOT CONVERGED!"); //strong convergence check
  }

  //obtain slip rate
  Real sliprate_mag = abs(solution);

  //update state variable
  // we apply separation of variables on first-order linear ODE of statevar
  Real coeff_LoverV = _rsf_L / sliprate_mag;
  Real coeff_exponent = exp(-sliprate_mag*_dt/_rsf_L);
  Real statevar_tplusdt = coeff_LoverV + (_statevar_old[_qp] - coeff_LoverV) * coeff_exponent;

  //*Compute shear traction at time t*
  Real T_mag = Tn * _rsf_a * asinh( 0.5*(sliprate_mag_old+sliprate_mag) * Z );

  ///Get Components
  Real T1 = T_mag * ( Ts_trial / Tmag_trial );
  Real T3 = T_mag * ( Td_trial / Tmag_trial );

  // Assign back traction in CZM
  RealVectorValue traction(T2 + _T2_o, -T1 + _T1_o, -T3 + _T3_o);
  _interface_traction[_qp] = traction;
  _dinterface_traction_djump[_qp] = 0;

  // Update statevar, slipratevar (strike, dip)
  _statevar[_qp] = statevar_tplusdt;
  _slipratevar[_qp](1) = sliprate_mag * ( Ts_trial / Tmag_trial );
  _slipratevar[_qp](2) = sliprate_mag * ( Td_trial / Tmag_trial );
}
