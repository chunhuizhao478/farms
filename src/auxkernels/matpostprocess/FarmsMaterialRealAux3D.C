//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FarmsMaterialRealAux3D.h"
#include "Assembly.h"

registerMooseObject("farmsApp", FarmsMaterialRealAux3D);

InputParameters
FarmsMaterialRealAux3D::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Creates a constant field in the domain.");
  params.addRequiredParam<std::string>("material_property_name", "The name of the material property to be retrieved: local_shear_jump, local_normal_jump, local_shear_traction, local_normal_traction, local_shear_jump_rate, local_normal_jump_rate, normal_x, normal_y, tangent_x, tangent_y");
  params.addRequiredCoupledVar("ini_shear_sts", "initial shear stress");
  params.addRequiredCoupledVar("ini_normal_sts", "initial normal stress");
  params.addParam<bool>("use_fractal_shear_stress", false, "use fractal shear stress");
  return params;
}

FarmsMaterialRealAux3D::FarmsMaterialRealAux3D(const InputParameters & parameters)
  : AuxKernel(parameters), 
  _material_property_name(getParam<std::string>("material_property_name")),
  _displacement_jump_global_x(getMaterialPropertyByName<Real>("jump_x")),
  _displacement_jump_global_y(getMaterialPropertyByName<Real>("jump_y")),
  _displacement_jump_global_z(getMaterialPropertyByName<Real>("jump_z")),
  _traction_global_x(getMaterialPropertyByName<Real>("traction_x")),
  _traction_global_y(getMaterialPropertyByName<Real>("traction_y")),
  _traction_global_z(getMaterialPropertyByName<Real>("traction_z")),
  _displacement_jump_global_x_old(getMaterialPropertyOldByName<Real>("jump_x")),
  _displacement_jump_global_y_old(getMaterialPropertyOldByName<Real>("jump_y")),
  _displacement_jump_global_z_old(getMaterialPropertyOldByName<Real>("jump_z")),
  _ini_shear_sts(coupledValue("ini_shear_sts")),
  _ini_normal_sts(coupledValue("ini_normal_sts")),
  _normals(_assembly.normals())
{
}

Real
FarmsMaterialRealAux3D::computeValue()
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
  Real jump_z = -1 * _displacement_jump_global_z[_qp];

  Real jump_x_old = -1 * _displacement_jump_global_x_old[_qp];
  Real jump_y_old = -1 * _displacement_jump_global_y_old[_qp];
  Real jump_z_old = -1 * _displacement_jump_global_z_old[_qp];
  
  /*
  * The traction has units: MPa (change, not total)
  */
  
  Real traction_x = -1 * _traction_global_x[_qp] / 1e6;
  Real traction_y = -1 * _traction_global_y[_qp] / 1e6;
  Real traction_z = -1 * _traction_global_z[_qp] / 1e6;

  /*
  * The normal is defined pointing at the positive side
  */

  Real normal_x = -1 * _normals[_qp](0);
  Real normal_y = -1 * _normals[_qp](1);
  Real normal_z = -1 * _normals[_qp](2);

  /*
  * The tangential direction is defined follows the direction of right-lateral motion
  */

  Real tangent_x = normal_y;
  Real tangent_y = -1 * normal_x;

  /*
  * We compute local jump, jump rate, traction in the local coordinate system
  */

  RealVectorValue global_jump(jump_x, jump_y, jump_z);
  RealVectorValue global_jump_old(jump_x_old, jump_y_old, jump_z_old);
  RealVectorValue global_traction(traction_x, traction_y, traction_z);
  RealVectorValue global_normal(normal_x, normal_y, normal_z);

  // Transform global jump and traction to local coordinate system
  Real local_jump_y = global_jump * global_normal;
  RealVectorValue local_jump_projectvec = global_jump - local_jump_y * global_normal;
  Real local_jump_x = local_jump_projectvec(0);
  Real local_jump_z = local_jump_projectvec(2);

  Real local_jump_y_old = global_jump_old * global_normal;
  RealVectorValue local_jump_projectvec_old = global_jump_old - local_jump_y_old * global_normal;
  Real local_jump_x_old = local_jump_projectvec_old(0);
  Real local_jump_z_old = local_jump_projectvec_old(2);

  Real local_traction_y = global_traction * global_normal;
  RealVectorValue local_traction_projectvec = global_traction - local_traction_y * global_normal;
  Real local_traction_x = local_traction_projectvec(0);
  Real local_traction_z = local_traction_projectvec(2);

  Real local_jump_rate_x = (local_jump_x - local_jump_x_old) / _dt;
  Real local_jump_rate_y = (local_jump_y - local_jump_y_old) / _dt;
  Real local_jump_rate_z = (local_jump_z - local_jump_z_old) / _dt;
  
  local_traction_x += _ini_shear_sts[_qp] / 1e6;
  local_traction_y += _ini_normal_sts[_qp] / 1e6;
  
  // Get the material property name and set the value accordingly
  if (_material_property_name == "local_shear_jump")
  {
    _value = local_jump_x;
  }
  else if (_material_property_name == "local_normal_jump")
  {
    _value = local_jump_y;
  }
  else if (_material_property_name == "local_dip_jump")
  {
    _value = local_jump_z;
  }
  else if (_material_property_name == "local_shear_jump_rate")
  {
    _value = local_jump_rate_x;
  }
  else if (_material_property_name == "local_normal_jump_rate")
  {
    _value = local_jump_rate_y;
  }
  else if (_material_property_name == "local_dip_jump_rate")
  {
    _value = local_jump_rate_z;
  }
  else if (_material_property_name == "local_shear_traction")
  {
    _value = local_traction_x;
  }
  else if (_material_property_name == "local_normal_traction")
  {
    _value = local_traction_y;
  }
  else if (_material_property_name == "local_dip_traction")
  {
    _value = local_traction_z;
  }
  else if (_material_property_name == "normal_x")
  {
    _value = normal_x;
  }
  else if (_material_property_name == "normal_y")
  {
    _value = normal_y;
  }
  else if (_material_property_name == "normal_z")
  {
    _value = normal_z;
  }
  else
  {
    mooseError("Invalid material property name: " + _material_property_name);
  }

  return _value;
}