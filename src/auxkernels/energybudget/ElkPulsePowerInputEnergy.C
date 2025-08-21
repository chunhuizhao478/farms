//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkPulsePowerInputEnergy.h"
#include "Assembly.h"
#include "Function.h"

registerMooseObject("farmsApp", ElkPulsePowerInputEnergy);

InputParameters
ElkPulsePowerInputEnergy::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Compute Pulse Power Input Energy");
  params.addParam<FunctionName>("function", -1, "Function being used to simulate pulse power");
  params.addParam<Real>("confinement_pressure", -1, "confinement pressure being used");
  params.addRequiredParam<int>("option","option: 1-function, 2-confinement_pressure");
  params.addRequiredCoupledVar("displacements", "The names of the displacement variables");
  params.addRequiredCoupledVar("old_displacements", "The names of the old displacement variables");
  return params;
}

ElkPulsePowerInputEnergy::ElkPulsePowerInputEnergy(const InputParameters & parameters)
  : AuxKernel(parameters),
    _option(getParam<int>("option")),
    _func(getFunction("function")),
    _confinement_pressure(getParam<Real>("confinement_pressure")),
    _ncomp(coupledComponents("displacements")),
    _normals(_assembly.normals())
{
  if (_ncomp != _mesh.dimension())
    paramError("displacements", "Number of entries must match the mesh dimension.");

  _disp.resize(_ncomp);
  _disp_old.resize(_ncomp);
  for (unsigned int j = 0; j < _ncomp; ++j){
    _disp[j] = &coupledValue("displacements", j); 
    _disp_old[j] = &coupledValue("old_displacements", j); 
  }

  if ( (_option != 1) and (_option != 2) ){
    mooseError("Must Specify Option in ElkPulsePowerInputEnergy!");
  }

}

Real
ElkPulsePowerInputEnergy::computeValue()
{
  
  //Evaluate normal displacements
  Real ddn = 0;
  for (unsigned int j = 0; j < _ncomp; ++j)
    ddn += _normals[_qp](j) * ( (*_disp[j])[_qp] - (*_disp_old[j])[_qp] );

  //return N/m based on option
  if (_option == 1){
    return (-_func.value(_t, _q_point[_qp])) * ddn;
  }
  else{
    return (-_confinement_pressure) * ddn;
  }

}