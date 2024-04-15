#include "InitialAlphaAD.h"

#include <random>

registerMooseObject("farmsApp", InitialAlphaAD);

InputParameters
InitialAlphaAD::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialAlphaAD::InitialAlphaAD(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialAlphaAD::value(Real /*t*/, const Point & p) const
{

  Real y_coord = p(1); //along the normal direction

  Real alpha_o = 0;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::weibull_distribution<double> wb_distribution(2.0,0.05);

  if (y_coord >= 0-1*0.1 and y_coord <= 0+1*0.1){
    alpha_o = 0.7;
  }
  else if (y_coord >= -45 and y_coord <= 45){
    alpha_o = wb_distribution(gen);
  }
  else{
    alpha_o = 0.0;
  }
  
  return alpha_o;

}