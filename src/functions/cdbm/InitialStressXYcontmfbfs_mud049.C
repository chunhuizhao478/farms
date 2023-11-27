#include "InitialStressXYcontmfbfs_mud049.h"

registerMooseObject("farmsApp", InitialStressXYcontmfbfs_mud049);

InputParameters
InitialStressXYcontmfbfs_mud049::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressXYcontmfbfs_mud049::InitialStressXYcontmfbfs_mud049(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressXYcontmfbfs_mud049::value(Real /*t*/, const Point & p) const
{

  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction

  Real T1_o = 0;

  if (x_coord >= -286 and x_coord <= 286 and y_coord >= -100 and y_coord <= 100)
  {
    //T1_o = 81.6e6;
    //raise overstress to 1%
    T1_o = 82.0524e6;
  }
  else{
    T1_o = 70e6;
  }
  
  return T1_o;

}