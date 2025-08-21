//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MonotoneIncreasingPostprocessor.h"

registerMooseObject("farmsApp", MonotoneIncreasingPostprocessor);

InputParameters
MonotoneIncreasingPostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredParam<PostprocessorName>("source", "The source postprocessor to monotone-ize");
  return params;
}

MonotoneIncreasingPostprocessor::MonotoneIncreasingPostprocessor(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _source(getPostprocessorValue("source")),
    _max_val(0.0),
    _initialized(false)
{
}

void MonotoneIncreasingPostprocessor::initialize()
{
  // no-op
}

void MonotoneIncreasingPostprocessor::execute()
{
  if (!_initialized)
  {
    _max_val = _source;
    _initialized = true;
  }
  else
    _max_val = std::max(_max_val, _source);
}

PostprocessorValue MonotoneIncreasingPostprocessor::getValue() const
{
  return _max_val;
}
