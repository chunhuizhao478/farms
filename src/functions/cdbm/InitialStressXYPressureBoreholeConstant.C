/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Will be Passed into AuxVariables
Chunhui Zhao
*/

#include "InitialStressXYPressureBoreholeConstant.h"

registerMooseObject("farmsApp", InitialStressXYPressureBoreholeConstant);

InputParameters
InitialStressXYPressureBoreholeConstant::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressXYPressureBoreholeConstant::InitialStressXYPressureBoreholeConstant(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressXYPressureBoreholeConstant::value(Real /*t*/, const Point & p) const
{
  
  //get coordinate of current point
  Real x_coord = p(0);
  Real y_coord = p(1);

  //compute ratio
  Real dx = x_coord - 0.0;
  Real dy = y_coord - 0.0;
  Real dl = std::sqrt(dx*dx+dy*dy);

  //compute pressure
  Real pval =0.0;
  if ( dl < 0.012 ){
    pval = 5e6;
  }

  return pval;

}