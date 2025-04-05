//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FarmsPorousFlowPermeabilityDamaged.h"

registerMooseObject("farmsApp", FarmsPorousFlowPermeabilityDamaged);

InputParameters
FarmsPorousFlowPermeabilityDamaged::validParams()
{
  InputParameters params = PorousFlowPermeabilityBase::validParams();
  params.addClassDescription(
      "This Material calculates the effective permeability tensor due to damage in smeared cracking model");
  return params;
}

FarmsPorousFlowPermeabilityDamaged::FarmsPorousFlowPermeabilityDamaged(const InputParameters & parameters)
  : PorousFlowPermeabilityBase(parameters),
    _effective_perm(getMaterialProperty<RealTensorValue>("effective_perm"))
{
}

void
FarmsPorousFlowPermeabilityDamaged::computeQpProperties()
{
  _permeability_qp[_qp] = _effective_perm[_qp];

  (*_dpermeability_qp_dvar)[_qp].assign(_num_var, RealTensorValue());
  (*_dpermeability_qp_dgradvar)[_qp].resize(LIBMESH_DIM);

  for (const auto i : make_range(Moose::dim))
    (*_dpermeability_qp_dgradvar)[_qp][i].assign(_num_var, RealTensorValue());

}