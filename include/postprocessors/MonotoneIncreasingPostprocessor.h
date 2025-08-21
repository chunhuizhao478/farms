//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralPostprocessor.h"

/**
 * Returns a non-decreasing running maximum of a source Postprocessor value.
 */
class MonotoneIncreasingPostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();
  MonotoneIncreasingPostprocessor(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual PostprocessorValue getValue() const override;

protected:
  const PostprocessorValue & _source;
  PostprocessorValue _max_val;
  bool _initialized;
};
