//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PoroCZMComputeGlobalTractionSmallStrain.h"

registerMooseObject("farmsApp", PoroCZMComputeGlobalTractionSmallStrain);

InputParameters
PoroCZMComputeGlobalTractionSmallStrain::validParams()
{
  InputParameters params = PoroCZMComputeGlobalTractionBase::validParams();

  params.addClassDescription(
      "Computes the czm traction in global coordinates for a small strain kinematic formulation");
  return params;
}

PoroCZMComputeGlobalTractionSmallStrain::PoroCZMComputeGlobalTractionSmallStrain(
    const InputParameters & parameters)
  : PoroCZMComputeGlobalTractionBase(parameters)
{
}

void
PoroCZMComputeGlobalTractionSmallStrain::computeEquilibriumTractionAndPressureAndDerivatives()
{
  _traction_global[_qp] = _czm_total_rotation[_qp] * _interface_traction[_qp];
  _dtraction_djump_global[_qp] = _czm_total_rotation[_qp] * _dinterface_traction_djump[_qp] *
                                 _czm_total_rotation[_qp].transpose();
  _dtraction_djump_global_vf[_qp] = _czm_total_rotation[_qp] * _dinterface_traction_djump_vf[_qp] *
                                 _czm_total_rotation[_qp].transpose();
  _dtraction_dpressure_global[_qp] = _czm_total_rotation[_qp] * _dinterface_traction_dpressure[_qp];                    
                                   
  _pressure_global[_qp] =  _interface_pressure[_qp];
  _dpressure_djump_global[_qp] = _czm_total_rotation[_qp] * _dinterface_pressure_djump[_qp];
  _dpressure_djump_global_vf[_qp] = _czm_total_rotation[_qp] * _dinterface_pressure_djump_vf[_qp];

                                                    
}