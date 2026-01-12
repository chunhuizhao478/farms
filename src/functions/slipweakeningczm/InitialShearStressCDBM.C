/*
Define Function for Initial Shear Stress for benchmark
*/

#include "InitialShearStressCDBM.h"
#include "SolutionUserObjectBase.h"
#include "MooseMesh.h"

#include <string.h>

registerMooseObject("farmsApp", InitialShearStressCDBM);

InputParameters
InitialShearStressCDBM::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>("peak_value", "Peak value of the initial shear stress");
  params.addRequiredParam<Real>("nucl_center_x", "X coordinate of the nucleation center");
  params.addRequiredParam<Real>("nucl_center_z", "Z coordinate of the nucleation center");
  params.addRequiredParam<Real>("nucl_size", "Size of the nucleation zone");
  params.addRequiredParam<Real>("elem_size", "Element size for the simulation");
  params.addRequiredParam<Real>("domain_value", "Domain value to use if not in nucleation zone");
  return params;
}

InitialShearStressCDBM::InitialShearStressCDBM(const InputParameters & parameters)
  : Function(parameters),
  _peak_value(getParam<Real>("peak_value")),
  _nucl_center_x(getParam<Real>("nucl_center_x")),
  _nucl_center_z(getParam<Real>("nucl_center_z")),
  _nucl_size(getParam<Real>("nucl_size")),
  _elem_size(getParam<Real>("elem_size")),
  _domain_value(getParam<Real>("domain_value"))
{
}

Real
InitialShearStressCDBM::value(Real t, const Point & p) const
{

  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the dip direction
  Real z_coord = p(2); //along the normal direction

  Real T1_o = 0.0;

  //get the solution object value
  Real domain_value = _domain_value;
    
  //tpv205
  if ((x_coord<=(_nucl_center_x+_nucl_size*0.5))&&(x_coord>=(_nucl_center_x-_nucl_size*0.5))&& (z_coord<=(_nucl_center_z+_nucl_size*0.5))&&(z_coord>=(_nucl_center_z-_nucl_size*0.5))&&(y_coord>=-1*_elem_size)&&(y_coord<=_elem_size))
  {
      T1_o = _peak_value;
  }
  else
  {
      T1_o = domain_value;
  }

  return T1_o;

}