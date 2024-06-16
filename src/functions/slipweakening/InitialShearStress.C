/*
Define Function for Initial Shear Stress for benchmark
*/

#include "InitialShearStress.h"

#include <string.h>

registerMooseObject("farmsApp", InitialShearStress);

InputParameters
InitialShearStress::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<std::string>("benchmark_type", "type of benchmark: tpv205, tpv14, tpv24");
  return params;
}

InitialShearStress::InitialShearStress(const InputParameters & parameters)
  : Function(parameters),
  _benchmark(getParam<std::string>("benchmark_type"))
{
}

Real
InitialShearStress::value(Real /*t*/, const Point & p) const
{

  Real x_coord = p(0); //along the x direction
  // Real y_coord = p(1); //along the y direction
  Real z_coord = p(2); //along the z direction

  Real T1_o = 0.0;

  //define option strings
  std::string tpv205 = "tpv205";
  std::string tpv14 = "tpv14";
  std::string tpv24 = "tpv24";

  if ( _benchmark == tpv205 ){
    
    //tpv205
    if ((x_coord<=(0.0+1.5e3))&&(x_coord>=(0.0-1.5e3))&& (z_coord<=(-7.5e3+1.5e3))&&(z_coord>=(-7.5e3-1.5e3)))
    {
        T1_o = 81.6e6;
    }
    else if ((x_coord<=(-7.5e3+1.5e3))&&(x_coord>=(-7.5e3-1.5e3)) && (z_coord<=(-7.5e3+1.5e3))&&(z_coord>=(-7.5e3-1.5e3)))
    {
        T1_o = 78.0e6;
    }
    else if ((x_coord<=(7.5e3+1.5e3))&&(x_coord>=(7.5e3-1.5e3)) && (z_coord<=(-7.5e3+1.5e3))&&(z_coord>=(-7.5e3-1.5e3)))
    {
        T1_o = 62.0e6;
    }
    else
    {
        T1_o = 70.0e6;
    }
  }
  else if ( _benchmark == tpv14 ){

    //tpv14
    if ((x_coord<=(-8e3+1.5e3))&&(x_coord>=(-8e3-1.5e3))&& (z_coord<=(-7.5e3+1.5e3))&&(z_coord>=(-7.5e3-1.5e3))){
        T1_o = 81.6e6;
    }
    else{
        T1_o = 70e6;
    }

  }
  else{
    mooseError("Unmatched benchmark type");
  }
  return T1_o;

}