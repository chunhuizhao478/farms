/*
AuxKernel for calculating total prinpal strain
*/

#include "PrincipalStrainCalc.h"

registerMooseObject("farmsApp", PrincipalStrainCalc);

InputParameters
PrincipalStrainCalc::validParams()
{
  InputParameters params = AuxKernel::validParams();
  
  params.addRequiredCoupledVar("eps_total_11", "total strain component 11"); 
  params.addRequiredCoupledVar("eps_total_12", "total strain component 12"); 
  params.addRequiredCoupledVar("eps_total_22", "total strain component 22"); 

  return params;
}

PrincipalStrainCalc::PrincipalStrainCalc(const InputParameters & parameters)
  : AuxKernel(parameters),
  _eps_total_11(coupledValue("eps_total_11")),
  _eps_total_12(coupledValue("eps_total_12")),
  _eps_total_22(coupledValue("eps_total_22")),
  _eps_total_11_old(coupledValueOld("eps_total_11")),
  _eps_total_12_old(coupledValueOld("eps_total_12")),
  _eps_total_22_old(coupledValueOld("eps_total_22"))
{
}

Real
PrincipalStrainCalc::computeValue()
{   
    //compute principal strain
    Real principal_strain = 0.5 * (_eps_total_11[_qp] + _eps_total_22[_qp]) + sqrt( pow( 0.5 * (_eps_total_11[_qp] - _eps_total_22[_qp]), 2) + pow( 0.5 * _eps_total_12[_qp], 2) );

    //compute old principal strain
    Real principal_strain_old = 0.5 * (_eps_total_11_old[_qp] + _eps_total_22_old[_qp]) + sqrt( pow( 0.5 * (_eps_total_11_old[_qp] - _eps_total_22_old[_qp]), 2) + pow( 0.5 * _eps_total_12_old[_qp], 2) );

    //compute principal strain rate
    Real principal_strain_rate = (principal_strain - principal_strain_old) / _dt;

    return principal_strain_rate;
}