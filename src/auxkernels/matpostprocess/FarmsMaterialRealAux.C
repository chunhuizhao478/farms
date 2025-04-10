//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FarmsMaterialRealAux.h"
#include "Assembly.h"

registerMooseObject("farmsApp", FarmsMaterialRealAux);

InputParameters
FarmsMaterialRealAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Creates a constant field in the domain.");
  params.addRequiredParam<std::string>("material_property_name", "The name of the material property to be retrieved: local_shear_jump, local_normal_jump, local_shear_traction, local_normal_traction, local_shear_jump_rate, local_normal_jump_rate, normal_x, normal_y, tangent_x, tangent_y");
  params.addRequiredCoupledVar("ini_shear_sts", "initial shear stress");
  params.addRequiredCoupledVar("ini_normal_sts", "initial normal stress");
  return params;
}

FarmsMaterialRealAux::FarmsMaterialRealAux(const InputParameters & parameters)
  : AuxKernel(parameters), 
  _material_property_name(getParam<std::string>("material_property_name")),
  _displacement_jump_global_x(getMaterialPropertyByName<Real>("jump_x")),
  _displacement_jump_global_y(getMaterialPropertyByName<Real>("jump_y")),
  _traction_global_x(getMaterialPropertyByName<Real>("traction_x")),
  _traction_global_y(getMaterialPropertyByName<Real>("traction_y")),
  _displacement_jump_global_x_old(getMaterialPropertyOldByName<Real>("jump_x")),
  _displacement_jump_global_y_old(getMaterialPropertyOldByName<Real>("jump_y")),
  _ini_shear_sts(coupledValue("ini_shear_sts")),
  _ini_normal_sts(coupledValue("ini_normal_sts")),
  _normals(_assembly.normals())
{
}

Real
FarmsMaterialRealAux::computeValue()
{

  Real _value = 0;

  // Get all the global variables
  
  /*
  * The CZM model defines the jump and traction as: negative side - positive side
  * The global variables are defined as: positive side - negative side
  * Here we reverse the sign to get the correct value
  */
  
  Real jump_x = -1 * _displacement_jump_global_x[_qp];
  Real jump_y = -1 * _displacement_jump_global_y[_qp];

  Real jump_x_old = -1 * _displacement_jump_global_x_old[_qp];
  Real jump_y_old = -1 * _displacement_jump_global_y_old[_qp];
  
  /*
  * The traction has units: MPa (change, not total)
  */
  
  Real traction_x = -1 * _traction_global_x[_qp] / 1e6;
  Real traction_y = -1 * _traction_global_y[_qp] / 1e6;

  /*
  * The normal is defined pointing at the positive side
  */

  Real normal_x = -1 * _normals[_qp](0);
  Real normal_y = -1 * _normals[_qp](1);

  /*
  * The tangential direction is defined follows the direction of right-lateral motion
  */

  Real tangent_x = normal_y;
  Real tangent_y = -1 * normal_x;

  /*
  * We compute local jump, jump rate, traction in the local coordinate system
  */

  Real local_jump_x = jump_x * tangent_x + jump_y * tangent_y;
  Real local_jump_y = jump_x * normal_x  + jump_y * normal_y;

  Real local_jump_x_old = jump_x_old * tangent_x + jump_y_old * tangent_y;
  Real local_jump_y_old = jump_x_old * normal_x  + jump_y_old * normal_y;

  Real local_jump_rate_x = (local_jump_x - local_jump_x_old) / _dt;
  Real local_jump_rate_y = (local_jump_y - local_jump_y_old) / _dt;

  Real local_traction_x = traction_x * tangent_x + traction_y * tangent_y;
  Real local_traction_y = traction_x * normal_x  + traction_y * normal_y;

  /*
  * now add initial stress state becomes total quantity
  */

  local_traction_x += _ini_shear_sts[_qp] / 1e6;
  local_traction_y += _ini_normal_sts[_qp] / 1e6;

  // 

  if (_material_property_name == "local_shear_jump")
  {
    _value = local_jump_x;
  }
  else if (_material_property_name == "local_normal_jump")
  {
    _value = local_jump_y;
  }
  else if (_material_property_name == "local_shear_jump_rate")
  {
    _value = local_jump_rate_x;
  }
  else if (_material_property_name == "local_normal_jump_rate")
  {
    _value = local_jump_rate_y;
  }
  else if (_material_property_name == "local_shear_traction")
  {
    _value = local_traction_x;
  }
  else if (_material_property_name == "local_normal_traction")
  {
    _value = local_traction_y;
  }
  else if (_material_property_name == "normal_x")
  {
    _value = normal_x;
  }
  else if (_material_property_name == "normal_y")
  {
    _value = normal_y;
  }
  else if (_material_property_name == "tangent_x")
  {
    _value = tangent_x;
  }
  else if (_material_property_name == "tangent_y")
  {
    _value = tangent_y;
  }
  else
  {
    mooseError("Invalid material property name: " + _material_property_name);
  }

  return _value;
}