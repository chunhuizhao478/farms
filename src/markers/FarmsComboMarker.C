//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FarmsComboMarker.h"
#include "FEProblemBase.h"

registerMooseObject("farmsApp", FarmsComboMarker);

InputParameters
FarmsComboMarker::validParams()
{
  InputParameters params = Marker::validParams();
  params.addRequiredParam<std::vector<MarkerName>>(
      "markers", "A list of marker names to combine into a single marker.");
  params.addRequiredParam<MarkerName>(
      "meshsize_marker", "The meshsize marker that acts as a gate - if not flagged, no refinement occurs.");
  params.addClassDescription("A marker that converts many markers into a single marker by "
                             "considering the maximum value of the listed markers (i.e., "
                             "refinement takes precedent). A special meshsize_marker acts as a gate - "
                             "if it's not flagged for refinement, no refinement occurs regardless of other markers.");
  return params;
}

FarmsComboMarker::FarmsComboMarker(const InputParameters & parameters)
  : Marker(parameters),
    _names(getParam<std::vector<MarkerName>>("markers")),
    _meshsize_marker_name(getParam<MarkerName>("meshsize_marker")),
    _block_restriction_mismatch(false),
    _meshsize_block_restriction_mismatch(false)
{
  // Set up regular markers
  for (const auto & marker_name : _names)
    _markers.push_back(&getMarkerValue(marker_name));

  // Set up meshsize marker
  _meshsize_marker = &getMarkerValue(_meshsize_marker_name);

  std::string other_block_restricted = "";
  for (const auto & marker_name : _names)
  {
    const auto var_ptr = &_subproblem.getVariable(_tid, marker_name);
    _marker_variables.push_back(var_ptr);

    // Check block restrictions
    if (blockIDs() != var_ptr->blockIDs())
      other_block_restricted += (other_block_restricted == "" ? "" : ", ") + marker_name;
  }

  // Check meshsize marker block restrictions
  _meshsize_marker_variable = &_subproblem.getVariable(_tid, _meshsize_marker_name);
  if (blockIDs() != _meshsize_marker_variable->blockIDs())
  {
    _meshsize_block_restriction_mismatch = true;
    other_block_restricted += (other_block_restricted == "" ? "" : ", ") + _meshsize_marker_name;
  }

  if (other_block_restricted != "")
  {
    _block_restriction_mismatch = true;
    paramInfo(
        "markers",
        "Combo marker and markers '" + other_block_restricted +
            "' do not share the same block restrictions. Markers outside their block restriction "
            "will not mark.");
  }
}

//enum MarkerValue
// {
//   DONT_MARK = -1,
//   COARSEN,      // 0
//   DO_NOTHING,   // 1
//   REFINE        // 2
// };

Marker::MarkerValue
FarmsComboMarker::computeElementMarker()
{
  // First check the meshsize marker - if it's not flagged for refinement, return DONT_MARK
  MarkerValue meshsize_value = DONT_MARK;
  
  if (!_meshsize_block_restriction_mismatch || 
      _meshsize_marker_variable->hasBlocks(_current_elem->subdomain_id()))
  {
    meshsize_value = static_cast<MarkerValue>((*_meshsize_marker)[0]);
  }
  
  // If meshsize marker is not flagged for refinement (DONT_MARK or COARSEN), don't refine
  if (meshsize_value == DO_NOTHING)
    return DONT_MARK;

  // Now process the regular markers - start with DONT_MARK because it's -1
  MarkerValue marker_value = DONT_MARK;

  // No need to check block restrictions if they all match
  if (!_block_restriction_mismatch)
    for (const auto & var : _markers)
      marker_value = std::max(marker_value, static_cast<MarkerValue>((*var)[0]));
  else
    for (const auto i : index_range(_markers))
      if (_marker_variables[i]->hasBlocks(_current_elem->subdomain_id()))
        marker_value = std::max(marker_value, static_cast<MarkerValue>((*_markers[i])[0]));

  // Return the combined marker value
  return marker_value;
}