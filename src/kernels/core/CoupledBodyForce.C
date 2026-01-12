//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CoupledBodyForce.h"

// MOOSE
#include "Function.h"
#include "Assembly.h"

registerMooseObject("farmsApp", CoupledBodyForce);
registerMooseObject("farmsApp", ADCoupledBodyForce);

template <bool is_ad>
InputParameters
CoupledBodyForceTempl<is_ad>::validParams()
{
  InputParameters params = GenericKernel<is_ad>::validParams();
  params.addClassDescription("Body force kernel scaled by a material property (e.g., density).");
  params.addParam<Real>("value", 1.0, "Coefficient to multiply by the body force term");
  params.addParam<FunctionName>("function", "1", "A function that describes the body force");
  params.addParam<PostprocessorName>(
      "postprocessor", 1, "A postprocessor whose value is multiplied by the body force");
  params.addRequiredParam<MaterialPropertyName>(
      "density_property_name", "Name of the material property to multiply (e.g., density)");
  params.declareControllable("value");
  return params;
}

template <bool is_ad>
CoupledBodyForceTempl<is_ad>::CoupledBodyForceTempl(const InputParameters & parameters)
  : GenericKernel<is_ad>(parameters),
    _scale(this->template getParam<Real>("value")),
    _function(getFunction("function")),
    _postprocessor(getPostprocessorValue("postprocessor")),
    _generic_q_point(this->_use_displaced_mesh ? &this->_assembly.template genericQPoints<is_ad>()
                                               : nullptr),
    _density(this->template getMaterialProperty<Real>(
        this->template getParam<MaterialPropertyName>("density_property_name")))
{
}

template <bool is_ad>
GenericReal<is_ad>
CoupledBodyForceTempl<is_ad>::computeQpResidual()
{
  const Real rho = _density[this->_qp];
  if (_generic_q_point)
    return -_test[_i][_qp] * _scale * _postprocessor * rho *
           _function.value(_t, (*_generic_q_point)[_qp]);
  else
    return -_test[_i][_qp] * _scale * _postprocessor * rho *
           _function.value(_t, _q_point[_qp]);
}

template class CoupledBodyForceTempl<false>;
template class CoupledBodyForceTempl<true>;
