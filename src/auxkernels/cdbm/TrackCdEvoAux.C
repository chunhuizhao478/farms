/*
AuxKernel of Passing Variable Time Derivative 
*/

#include "TrackCdEvoAux.h"

registerMooseObject("farmsApp", TrackCdEvoAux);

InputParameters
TrackCdEvoAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  
  params.addRequiredParam<Real>(   "C_d_min", "coefficient gives positive damage evolution (small strain e < 1e-4 threshold value)");
  params.addRequiredParam<Real>("mechanical_strain_rate_threshold", "threshold value for strain rate such that Cd takes constant value Cd_min if strain rate below this value.");
  params.addRequiredCoupledVar("mechanical_strain_rate", "strain rate"); 
  params.addRequiredParam<Real>(     "m", "Cd power-law correction index");

  return params;
}

TrackCdEvoAux::TrackCdEvoAux(const InputParameters & parameters)
  : AuxKernel(parameters),
  _Cd_min(getParam<Real>("C_d_min")),
  _mechanical_strain_rate_threshold(getParam<Real>("mechanical_strain_rate_threshold")),
  _mechanical_strain_rate(coupledValue("mechanical_strain_rate")),
  _m(getParam<Real>("m"))
{
}

Real
TrackCdEvoAux::computeValue()
{
    //Power-law correction
    //Initialize Cd
    Real Cd = 0;
    //power-law correction on coefficient Cd(function of strain rate)
    if ( _mechanical_strain_rate[_qp] < _mechanical_strain_rate_threshold ) //Cd remain constant
    {
        Cd = _Cd_min;
    }
    else{ //Cd follows power-law
        Cd = pow(10, 1 + _m * log10( _mechanical_strain_rate[_qp] / _mechanical_strain_rate_threshold ) ) * _Cd_min;
    }
    return Cd;
}