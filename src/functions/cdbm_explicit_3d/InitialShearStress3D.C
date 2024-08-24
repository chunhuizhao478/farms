#include "InitialShearStress3D.h"

registerMooseObject("farmsApp", InitialShearStress3D);

InputParameters
InitialShearStress3D::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<std::vector<Real>>("nucl_center", "nucleation center (x,y,z)");
  params.addRequiredParam<Real>("e_sigma","the standard deviation used in apply damage value normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("min_val","min value");
  params.addRequiredParam<Real>("max_val","max value");
  return params;
}

InitialShearStress3D::InitialShearStress3D(const InputParameters & parameters)
  : Function(parameters),
  _nucl_center(getParam<std::vector<Real>>("nucl_center")),
  _sigma(getParam<Real>("e_sigma")),
  _min_val(getParam<Real>("min_val")),
  _max_val(getParam<Real>("max_val"))
{
}

Real
InitialShearStress3D::value(Real /*t*/, const Point & p) const
{
  Real xcoord = p(0); //along the strike direction
  Real ycoord = p(1); //along the dip direction
  Real zcoord = p(2); //along the normal direction

  Real T1_o = 0;

  Real r = std::sqrt(pow(xcoord - _nucl_center[0],2) + pow(ycoord - _nucl_center[1],2) + pow(zcoord - _nucl_center[2],2));
  T1_o = _min_val + std::max((_max_val-_min_val) * exp(-1.0*(std::pow(r,2))/(_sigma*_sigma)),0.0);
  
  return T1_o;

}