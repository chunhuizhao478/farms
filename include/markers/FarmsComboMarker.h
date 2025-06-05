//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Marker.h"

/**
 * Combines multiple marker fields.  The most conservative (requesting additional refinement) wins.
 * A special meshsize_marker acts as a gate - if it's not flagged, no refinement occurs.
 */
class FarmsComboMarker : public Marker
{
public:
  static InputParameters validParams();

  FarmsComboMarker(const InputParameters & parameters);

protected:
  virtual MarkerValue computeElementMarker() override;

  /// Names of the markers contributing to the combo
  const std::vector<MarkerName> & _names;

  /// Name of the meshsize marker that acts as a gate
  const MarkerName & _meshsize_marker_name;

  /// Pointers to the markers contributing to the Combo marker
  std::vector<const VariableValue *> _markers;

  /// Pointer to the meshsize marker
  const VariableValue * _meshsize_marker;

  /// Boolean to keep track of whether any marker does not have the same block restriction
  bool _block_restriction_mismatch;

  /// Boolean to keep track of whether meshsize marker has different block restriction
  bool _meshsize_block_restriction_mismatch;

  /// Pointers to the variables for the markers
  std::vector<const MooseVariableFieldBase *> _marker_variables;

  /// Pointer to the variable for the meshsize marker
  const MooseVariableFieldBase * _meshsize_marker_variable;
};