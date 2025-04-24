#pragma once

#include "AuxKernel.h"

/**
 * Copies one variable onto an auxiliary variable.
 * Now stores coupled data as pointers to allow proper construction.
 */
class CheckAlpha : public AuxKernel
{
public:
  static InputParameters validParams();

  CheckAlpha(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  // Instead of reference members, use pointers so that construction and assignment work.
  const VariableValue * _v;
  const MooseVariable * _source_variable;
  const VariableValue * _initial_damage_aux;

  // We still store the state as before.
  unsigned short _state;
};