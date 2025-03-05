/*
Define Function for Initial Shear Stress for benchmark
*/

#include "InitialShearStressCDBM.h"

#include <string.h>

registerMooseObject("farmsApp", InitialShearStressCDBM);

InputParameters
InitialShearStressCDBM::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<std::string>("benchmark_type", "type of benchmark: tpv205, tpv14, tpv24");
  return params;
}

InitialShearStressCDBM::InitialShearStressCDBM(const InputParameters & parameters)
  : Function(parameters),
  _benchmark(getParam<std::string>("benchmark_type"))
{
}

Real
InitialShearStressCDBM::value(Real /*t*/, const Point & p) const
{

  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the dip direction
  Real z_coord = p(2); //along the normal direction

  Real T1_o = 0.0;

  //define option strings
  std::string tpv205 = "tpv205";
  std::string tpv14 = "tpv14";
  std::string tpv24 = "tpv24";

  if ( _benchmark == tpv205 ){
    
    //tpv205
    if ((x_coord<=(0.0+1.5e3))&&(x_coord>=(0.0-1.5e3))&& (y_coord<=(-7.5e3+1.5e3))&&(y_coord>=(-7.5e3-1.5e3))&&(z_coord>=-200)&&(z_coord<=200))
    {
        T1_o = 81.6e6;
    }
    else
    {
        T1_o = 60.0e6;
    }
  }
  else if ( _benchmark == tpv24 ){

    //tpv14
    if ((x_coord<=(0+1.5e3))&&(x_coord>=(0-1.5e3))&& (y_coord<=(-8e3+1.5e3))&&(y_coord>=(-8e3-1.5e3))&&(z_coord<=(0+200))&&(z_coord>=(0-200))){
        T1_o = 60.58e6;
    }
    else{
        if(-y_coord<15600){
          T1_o = 1 * (-0.169029 * ( (-2700 * 9.8 * (-y_coord)) + (1000 * 9.8 * (-y_coord)) ));
        }
        else{
          T1_o = 0.0;
        }
    }

  }
  else{
    mooseError("Unmatched benchmark type");
  }
  return T1_o;

}