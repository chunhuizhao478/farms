//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ConfigurationalForceAux.h"
#include "ConfigurationalForceUserObject.h"

registerMooseObject("farmsApp", ConfigurationalForceAux);

InputParameters
ConfigurationalForceAux::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addClassDescription(
      "Outputs a component of the configurational force vector computed by "
      "ConfigurationalForceUserObject to an AuxVariable for visualization in Paraview/Exodus.");

  params.addRequiredParam<UserObjectName>(
      "configurational_force_uo",
      "Name of the ConfigurationalForceUserObject that computed the forces");

  MooseEnum component("x=0 y=1 z=2 magnitude=3", "magnitude");
  params.addParam<MooseEnum>("component",
                             component,
                             "Which component to output: x, y, z, or magnitude (default)");

  return params;
}

ConfigurationalForceAux::ConfigurationalForceAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _config_force_uo(getUserObject<ConfigurationalForceUserObject>("configurational_force_uo")),
    _component(getParam<MooseEnum>("component"))
{
}

Real
ConfigurationalForceAux::computeValue()
{
  // Get the node ID of the current node
  dof_id_type node_id = _current_node->id();

  // Query the UserObject for the appropriate component or magnitude
  if (_component == 3) // magnitude
    return _config_force_uo.getForceMagnitude(node_id);
  else // x, y, or z component
    return _config_force_uo.getForceComponent(node_id, _component);
}
