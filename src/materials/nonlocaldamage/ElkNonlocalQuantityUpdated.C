//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkNonlocalQuantityUpdated.h"
#include "FEProblem.h"

registerMooseObject("farmsApp", ElkNonlocalQuantityUpdated);

InputParameters
ElkNonlocalQuantityUpdated::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Retrieve nonlocal averaged property from radial averaging user object. "
      "This generic material can handle any scalar material property that needs nonlocal averaging "
      "(e.g., strain invariant ratio 'xi', deviatoric strain rate, etc.)");

  params.addRequiredParam<UserObjectName>("average_UO",
                                          "Radial Average user object (ElkRadialAverageUpdated)");

  // Make property names configurable for generality
  params.addParam<MaterialPropertyName>("nonlocal_property_name",
                                        "eqstrain_nonlocal",
                                        "Name of the nonlocal averaged property to declare "
                                        "(default: eqstrain_nonlocal for backward compatibility)");

  params.addParam<MaterialPropertyName>("initial_property_name",
                                        "eqstrain_nonlocal_initial",
                                        "Name of the initial value property to retrieve "
                                        "(default: eqstrain_nonlocal_initial for backward compatibility)");

  return params;
}

ElkNonlocalQuantityUpdated::ElkNonlocalQuantityUpdated(const InputParameters & parameters)
  : Material(parameters),
    _prop_nonlocal(declareProperty<Real>(getParam<MaterialPropertyName>("nonlocal_property_name"))),
    _prop_nonlocal_old(getMaterialPropertyOld<Real>(getParam<MaterialPropertyName>("nonlocal_property_name"))),
    _average(getUserObject<ElkRadialAverageUpdated>("average_UO").getAverage()),
    _prev_elem(nullptr),
    _prop_initial(getMaterialProperty<Real>(getParam<MaterialPropertyName>("initial_property_name"))),
    _step(_fe_problem.timeStep())
{
}

//Rules: See https://github.com/idaholab/moose/discussions/19450
//Only the object that declares the material property can assign values to it.
//Objects can request material properties, gaining read-only access to their values.
//When any object (including the object that declares it) requests the old value of a material property, that property becomes "stateful".
//All stateful material properties must be initialized within the initQpStatefulProperties call.
//
void
ElkNonlocalQuantityUpdated::initQpStatefulProperties()
{
  _prop_nonlocal[_qp] = _prop_initial[_qp];
}

void
ElkNonlocalQuantityUpdated::computeQpProperties()
{
  // For the first time step, use initial value
  if (_step == 1)
  {
    _prop_nonlocal[_qp] = _prop_initial[_qp];
    return;
  }

  // Now update the nonlocal property from the averaging user object
  // Only update iterator when we change to another element. This is for
  // computational efficiency related to map lookup.
  if (_prev_elem != _current_elem)
  {
    _average_iterator = _average.find(_current_elem->id());
    _prev_elem = _current_elem;
  }

  // Check that we found the element
  if (_average_iterator != _average.end())
    _prop_nonlocal[_qp] = _average_iterator->second[_qp];
  else
    // During startup the map is not made yet or
    // if AMR is used then the new element will not be found but it should
    // already have an old nonlocal value that needs to be preserved
    _prop_nonlocal[_qp] = _prop_nonlocal_old[_qp];
}
