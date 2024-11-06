//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StaticStressCycleSim.h"

/**
 *  Created by Chunhui Zhao, Nov 5th, 2024
 *  Material used in extracting the static stress after the first solve
 */
registerMooseObject("farmsApp", StaticStressCycleSim);

InputParameters
StaticStressCycleSim::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in defining initial damage profile");
  return params;
}

StaticStressCycleSim::StaticStressCycleSim(const InputParameters & parameters)
  : Material(parameters),
  _static_stress(declareProperty<RankTwoTensor>("static_stress")),
  _static_stress_old(getMaterialPropertyOldByName<RankTwoTensor>("static_stress")),
  _pk2_stress(getMaterialProperty<RankTwoTensor>("pk2_stress"))
{
}

void
StaticStressCycleSim::initQpStatefulProperties()
{
  _static_stress[_qp].zero();
}

void
StaticStressCycleSim::computeQpProperties()
{
  // Debug output
  // mooseInfo("Current PK2 stress: ", _pk2_stress[_qp]);
  // mooseInfo("Current static stress: ", _static_stress[_qp]);

  // Store PK2 stress after first timestep is solved
  if (_t_step == 1 || _t_step == 2)
  {
    _static_stress[_qp] = _pk2_stress[_qp];
    //mooseInfo("Storing PK2 stress after first step solved, time = ", _t);
  }
  else
  {
    _static_stress[_qp] = _static_stress_old[_qp];
  }

  //mooseInfo("Final static stress: ", _static_stress[_qp]);  
}