#include "InitialStressXYcontmfbfs.h"

registerMooseObject("farmsApp", InitialStressXYcontmfbfs);

InputParameters
InitialStressXYcontmfbfs::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressXYcontmfbfs::InitialStressXYcontmfbfs(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressXYcontmfbfs::value(Real /*t*/, const Point & p) const
{

  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction

  Real T1_o = 0;

  if (x_coord >= -100 and x_coord <= 100 and y_coord >= -100 and y_coord <= 100)
  {
    T1_o = 81.6e6;
  }
  else{
    T1_o = 70e6;
  }
  
  return T1_o;

}