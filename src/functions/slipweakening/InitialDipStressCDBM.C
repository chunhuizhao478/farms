/*
Define Function for Initial Shear Stress for benchmark
*/

#include "InitialDipStressCDBM.h"

#include <string.h>

registerMooseObject("farmsApp", InitialDipStressCDBM);

InputParameters
InitialDipStressCDBM::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<std::string>("benchmark_type", "type of benchmark: tpv205, tpv14, tpv24");
  return params;
}

InitialDipStressCDBM::InitialDipStressCDBM(const InputParameters & parameters)
  : Function(parameters),
  _benchmark(getParam<std::string>("benchmark_type"))
{
}

Real
InitialDipStressCDBM::value(Real /*t*/, const Point & p) const
{

  //Real x_coord = p(0); //along the x direction
  Real y_coord = p(1); //along the y direction
  //Real z_coord = p(2); //along the z direction

  Real T1_o = 0.0;

  //define option strings
  std::string tpv205 = "tpv205";
  std::string tpv14 = "tpv14";
  std::string tpv24 = "tpv24";

  if ( _benchmark == tpv205 ){
    
    //tpv205
    T1_o = -2670 * 9.8 * abs(y_coord);
  }
  else if ( _benchmark == tpv14 ){
  }
  else{
    mooseError("Unmatched benchmark type");
  }
  return T1_o;

}