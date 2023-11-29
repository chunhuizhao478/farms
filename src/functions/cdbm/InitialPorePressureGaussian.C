#include "InitialPorePressureGaussian.h"

registerMooseObject("farmsApp", InitialPorePressureGaussian);

InputParameters
InitialPorePressureGaussian::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialPorePressureGaussian::InitialPorePressureGaussian(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialPorePressureGaussian::value(Real /*t*/, const Point & p) const
{

  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction

  Real T1_o = 0;

  if (x_coord >= -200 and x_coord <= 200 and y_coord >= -100 and y_coord <= 100)
  {
    //T1_o = 81.6e6;
    //raise overstress to 1%
    T1_o = 70e6;
  }
  else{
    T1_o = 70e6;
  }
  
  return T1_o;

}