/*
AuxKernel of Passing Variable
*/

#include "ComputePlasticWorkRate.h"

registerMooseObject("farmsApp", ComputePlasticWorkRate);

InputParameters
ComputePlasticWorkRate::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("eps_p_11_aux","plastic strain 11");
  params.addRequiredCoupledVar("eps_p_12_aux","plastic strain 12");
  params.addRequiredCoupledVar("eps_p_22_aux","plastic strain 22");
  params.addRequiredCoupledVar("eps_p_33_aux","plastic strain 33");
  params.addRequiredCoupledVar("sts_change_11_aux","stress change 11");
  params.addRequiredCoupledVar("sts_change_12_aux","stress change 12");
  params.addRequiredCoupledVar("sts_change_22_aux","stress change 22");
  params.addRequiredCoupledVar("sts_change_33_aux","stress change 33");
  params.addRequiredCoupledVar("sts_initial_11_aux","initial stress 11");
  params.addRequiredCoupledVar("sts_initial_12_aux","initial stress 12");
  params.addRequiredCoupledVar("sts_initial_22_aux","initial stress 22");
  params.addRequiredCoupledVar("sts_initial_33_aux","initial stress 33");
  return params;
}

ComputePlasticWorkRate::ComputePlasticWorkRate(const InputParameters & parameters)
  : AuxKernel(parameters),
  
  _eps_p_11(coupledValue("eps_p_11_aux")),
  _eps_p_11_old(coupledValueOld("eps_p_11_aux")),
  _eps_p_12(coupledValue("eps_p_12_aux")),
  _eps_p_12_old(coupledValueOld("eps_p_12_aux")),
  _eps_p_22(coupledValue("eps_p_22_aux")),
  _eps_p_22_old(coupledValueOld("eps_p_22_aux")),
  _eps_p_33(coupledValue("eps_p_33_aux")),
  _eps_p_33_old(coupledValueOld("eps_p_33_aux")),
  _sts_change_11(coupledValue("sts_change_11_aux")),
  _sts_change_12(coupledValue("sts_change_12_aux")),
  _sts_change_22(coupledValue("sts_change_22_aux")),
  _sts_change_33(coupledValue("sts_change_33_aux")),
  _sts_initial_11(coupledValue("sts_initial_11_aux")),
  _sts_initial_12(coupledValue("sts_initial_12_aux")),
  _sts_initial_22(coupledValue("sts_initial_22_aux")),
  _sts_initial_33(coupledValue("sts_initial_33_aux"))

{
}

Real
ComputePlasticWorkRate::computeValue()
{

  //Compute Inelastic Strain Rate
  Real eps_p_11_rate = ( _eps_p_11[_qp] - _eps_p_11_old[_qp] ) / _dt;
  Real eps_p_12_rate = ( _eps_p_12[_qp] - _eps_p_12_old[_qp] ) / _dt;
  Real eps_p_22_rate = ( _eps_p_22[_qp] - _eps_p_22_old[_qp] ) / _dt;
  Real eps_p_33_rate = ( _eps_p_33[_qp] - _eps_p_33_old[_qp] ) / _dt;

  //Compute total stress
  Real sts_total_11 = _sts_change_11[_qp] + _sts_initial_11[_qp];
  Real sts_total_12 = _sts_change_12[_qp] + _sts_initial_12[_qp];
  Real sts_total_22 = _sts_change_22[_qp] + _sts_initial_22[_qp];
  Real sts_total_33 = _sts_change_33[_qp] + _sts_initial_33[_qp];

  //Compute Plastic Work Rate
  Real plastic_work_rate = eps_p_11_rate * sts_total_11 + 2 * eps_p_12_rate * sts_total_12 + eps_p_22_rate * sts_total_22 + eps_p_33_rate * sts_total_33;

  return plastic_work_rate;
}