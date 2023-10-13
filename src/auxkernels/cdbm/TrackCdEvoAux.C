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
  params.addParam<Real>( "scale", 1.0, "scale the Cd power-law");
  params.addRequiredParam<int>( "option", "option 1 : Cd power-law; option 2 : use constant Cd");
  params.addParam<Real>( "Cd_constant", 0.0, "constant Cd value for option 2 only");

  return params;
}

TrackCdEvoAux::TrackCdEvoAux(const InputParameters & parameters)
  : AuxKernel(parameters),
  _Cd_min(getParam<Real>("C_d_min")),
  _mechanical_strain_rate_threshold(getParam<Real>("mechanical_strain_rate_threshold")),
  _mechanical_strain_rate(coupledValue("mechanical_strain_rate")),
  _m(getParam<Real>("m")),
  _scale(getParam<Real>("scale")),
  _option(getParam<int>("option")),
  _Cd_constant(getParam<Real>("Cd_constant"))
{
}

Real
TrackCdEvoAux::computeValue()
{   
    Real Cd = 0.0;
    Real Cd_bound = 400;
    //Check options
    if ( _option == 1 ){

      //power-law correction on coefficient Cd(function of strain rate)
      if ( _mechanical_strain_rate[_qp] > _mechanical_strain_rate_threshold ) //Cd follows power-law
      { 
        if ( Cd < Cd_bound ){
          Cd = std::min(_scale * pow(10, 1 + _m * log10( _mechanical_strain_rate[_qp] / _mechanical_strain_rate_threshold ) ) * _Cd_min, Cd_bound);
        }
        else{
          Cd = Cd_bound;
        }
      }
      else if ( _mechanical_strain_rate[_qp] > 0 ){ //Cd remains constant
        Cd = _Cd_min;
      }

    }
    else if ( _option == 2 ){

      if ( _Cd_constant == 0.0 ){
        mooseError("For option 2, need to provide nonzero Cd_constant value !");
      }
      else{
        Cd = _Cd_constant;
      }

    }
    else{
      mooseError("Please provide valid option number!");
    }
    return Cd;
}