#include "InitialNormalStressAD.h"

#include <random>

registerMooseObject("farmsApp", InitialNormalStressAD);

InputParameters
InitialNormalStressAD::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialNormalStressAD::InitialNormalStressAD(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialNormalStressAD::value(Real /*t*/, const Point & p) const
{

  Real x_coord = p(0); //along the normal direction
  Real y_coord = p(1); //along the normal direction

  // Real T2_o_perturb = 1e5;

  std::random_device rd;
  std::mt19937 gen(rd());
  // std::weibull_distribution<double> wb_distribution(2.0,T2_o_perturb);
  std::uniform_real_distribution<> dis(1e4, 1e5);

  Real T2_o = 0.0;

  if (y_coord >= 0-1*0.1 and y_coord <= 0+1*0.1){
    if (x_coord >= 0-1*0.05 and x_coord <= 0+1*0.05 and y_coord >= 0-1*0.05 and y_coord <= 0+1*0.05){
      T2_o = -50e6 + 0.5e6;
    }
  }
  else{
    T2_o = 0.0;
  }
  
  return T2_o;

}