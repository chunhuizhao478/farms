//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkNonlocalEqstrain.h"

registerMooseObject("farmsApp", ElkNonlocalEqstrain);

InputParameters
ElkNonlocalEqstrain::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Compute the local elastic energy");
  params.addRequiredParam<UserObjectName>("average_UO", "Radial Average user object");
  return params;
}

ElkNonlocalEqstrain::ElkNonlocalEqstrain(const InputParameters & parameters)
  : Material(parameters),
    _eqstrain_nonlocal(declareProperty<Real>("eqstrain_nonlocal")),
    _eqstrain_nonlocal_old(getMaterialPropertyOld<Real>("eqstrain_nonlocal")),
    _average(getUserObject<ElkRadialAverage>("average_UO").getAverage()),
    _prev_elem(nullptr)
{
}

//Rules:See https://github.com/idaholab/moose/discussions/19450
//Only the object that declares the material property can assign values to it.
//Objects can request material properties, gaining read-only access to their values.
//When any object (including the object that declares it) requests the old value of a material property, that property becomes "stateful".
//All stateful material properties must be initialized within the initQpStatefulProperties call. 
//
void
ElkNonlocalEqstrain::initQpStatefulProperties()
{
  _eqstrain_nonlocal[_qp] = 0.0;
}

void
ElkNonlocalEqstrain::computeQpProperties()
{

  // Now update the nonlocal damage model
  _average_eqstrain_nonlocal = _average.find(_current_elem->id());
  _eqstrain_nonlocal[_qp] = _average_eqstrain_nonlocal->second[_qp];

  // Now update the nonlocal damage model
  //Only update iterator when we change to another element. This is for
  //computational costs related to map lookup.
  // if (_prev_elem != _current_elem)
  // {
  //   _average_eqstrain_nonlocal = _average.find(_current_elem->id());
  //   _prev_elem = _current_elem;
  // }
  // // Check that we found the new element
  // if (_average_eqstrain_nonlocal != _average.end())
  //   _eqstrain_nonlocal[_qp] = _average_eqstrain_nonlocal->second[_qp];
  // else
  //   // during startup the map is not made yet or
  //   // if AMR is used then the new element will not be found but it should
  //   // already have an old nonlocal damage value that needs to perserved
  //   _eqstrain_nonlocal[_qp] = _eqstrain_nonlocal_old[_qp];

}
