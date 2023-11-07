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
  params.addRequiredCoupledVar("B","breakage parameter");
  params.addRequiredParam<Real>(             "C_g", "material parameter: compliance or fluidity of the fine grain granular material");
  params.addRequiredParam<Real>(              "m1", "coefficient of power law indexes");
  params.addRequiredParam<Real>(              "m2", "coefficient of power law indexes");
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
  _sts_initial_33(coupledValue("sts_initial_33_aux")),
  _B(coupledValue("B")),
  _C_g(getParam<Real>("C_g")),
  _m1(getParam<Real>("m1")),
  _m2(getParam<Real>("m2"))
{
}

Real
ComputePlasticWorkRate::computeValue()
{

  //Compute total stress
  Real sts_total_11 = _sts_change_11[_qp] + _sts_initial_11[_qp];
  Real sts_total_12 = _sts_change_12[_qp] + _sts_initial_12[_qp];
  Real sts_total_22 = _sts_change_22[_qp] + _sts_initial_22[_qp];
  Real sts_total_33 = _sts_change_33[_qp] + _sts_initial_33[_qp];

  //Compute deviatroic stress components
  Real sts_total_tr = sts_total_11 + sts_total_22 + sts_total_33;
  Real sts_d_11 = sts_total_11 - 1/3 * sts_total_tr;
  Real sts_d_12 = sts_total_12;
  Real sts_d_22 = sts_total_22 - 1/3 * sts_total_tr;
  Real sts_d_33 = sts_total_33 - 1/3 * sts_total_tr;
 
  //Compute Inelastic Strain Rate
  Real eps_p_11_rate = _C_g * pow( _B[_qp], _m1 ) * pow( sts_d_11 , _m2 );
  Real eps_p_12_rate = _C_g * pow( _B[_qp], _m1 ) * pow( sts_d_12 , _m2 );
  Real eps_p_22_rate = _C_g * pow( _B[_qp], _m1 ) * pow( sts_d_22 , _m2 );
  Real eps_p_33_rate = _C_g * pow( _B[_qp], _m1 ) * pow( sts_d_33 , _m2 );

  //Compute Plastic Work Rate
  Real plastic_work_rate = eps_p_11_rate * sts_total_11 + 2 * eps_p_12_rate * sts_total_12 + eps_p_22_rate * sts_total_22 + eps_p_33_rate * sts_total_33;

  return plastic_work_rate;
}