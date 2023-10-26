#include "InitialStressXYcontmfbfs_mud04.h"

registerMooseObject("farmsApp", InitialStressXYcontmfbfs_mud04);

InputParameters
InitialStressXYcontmfbfs_mud04::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressXYcontmfbfs_mud04::InitialStressXYcontmfbfs_mud04(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressXYcontmfbfs_mud04::value(Real /*t*/, const Point & p) const
{

  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction

  Real T1_o = 0;

  if (x_coord >= -200 and x_coord <= 200 and y_coord >= -100 and y_coord <= 100)
  {
    T1_o = 81.6e6;
  }
  else{
    T1_o = 70e6;
  }
  
  return T1_o;

}