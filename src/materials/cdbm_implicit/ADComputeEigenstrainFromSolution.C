//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeEigenstrainFromSolution.h"
#include "RankTwoTensor.h"
#include "Function.h"
#include "Conversion.h" // for stringify

registerMooseObject("farmsApp", ADComputeEigenstrainFromSolution);

InputParameters
ADComputeEigenstrainFromSolution::validParams()
{
  InputParameters params = ADComputeEigenstrainBase::validParams();
  params.addClassDescription("Computes an eigenstrain from an initial strain");
  params.addRequiredParam<std::vector<FunctionName>>(
      "initial_strain",
      "A list of functions describing the initial strain.  There must be 9 of these, corresponding "
      "to the xx, yx, zx, xy, yy, zy, xz, yz, zz components respectively.  To compute the "
      "eigenstrain correctly, your elasticity tensor should not be time-varying in the first "
      "timestep");
  params.addCoupledVar("initial_strain_aux",
                       "A list of 9 AuxVariables describing the initial strain.  If provided, each "
                       "of these is multiplied by its corresponding initial_strain function to "
                       "obtain the relevant component of initial strain.");
  params.addParam<std::string>("base_name",
                               "The base_name for the elasticity tensor that will be "
                               "used to compute strain from strain.  Do not provide "
                               "any base_name if your elasticity tensor does not use "
                               "one.");
  return params;
}

ADComputeEigenstrainFromSolution::ADComputeEigenstrainFromSolution(
    const InputParameters & parameters)
  : ADComputeEigenstrainBase(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _eigenstrain_old(getMaterialPropertyOld<RankTwoTensor>(_eigenstrain_name)),
    _ini_aux_provided(isParamValid("initial_strain_aux")),
    _ini_aux(_ini_aux_provided ? coupledValues("initial_strain_aux")
                               : std::vector<const VariableValue *>())
{
  const std::vector<FunctionName> & fcn_names(
      getParam<std::vector<FunctionName>>("initial_strain"));
  const std::size_t num = fcn_names.size();

  if (num != LIBMESH_DIM * LIBMESH_DIM)
    paramError(
        "initial_strain",
        "ADComputeEigenstrainFromSolution: " + Moose::stringify(LIBMESH_DIM * LIBMESH_DIM) +
            " initial strain functions must be provided.  You supplied " + Moose::stringify(num) +
            "\n");

  _initial_strain_fcn.resize(num);
  for (unsigned i = 0; i < num; ++i)
    _initial_strain_fcn[i] = &getFunctionByName(fcn_names[i]);

  if (_ini_aux_provided)
  {
    const std::size_t aux_size = coupledComponents("initial_strain_aux");
    if (aux_size != LIBMESH_DIM * LIBMESH_DIM)
      paramError("initial_strain_aux",
                 "ADComputeEigenstrainFromSolution: If you supply initial_strain_aux, " +
                     Moose::stringify(LIBMESH_DIM * LIBMESH_DIM) +
                     " values must be given.  You supplied " + Moose::stringify(aux_size) + "\n");
  }
}

void
ADComputeEigenstrainFromSolution::computeQpEigenstrain()
{
  if (_t_step == 0) //load the solution from initial
  {
    ADRankTwoTensor initial_strain;
    for (unsigned i = 0; i < LIBMESH_DIM; ++i)
      for (unsigned j = 0; j < LIBMESH_DIM; ++j)
      {
        initial_strain(i, j) = _initial_strain_fcn[i * LIBMESH_DIM + j]->value(_t, _q_point[_qp]);
        if (_ini_aux_provided)
          initial_strain(i, j) *= (*_ini_aux[i * LIBMESH_DIM + j])[_qp];
      }
    
    //invSymm only works for non-AD
    _eigenstrain[_qp] = initial_strain;
  }
  else
    _eigenstrain[_qp] = _eigenstrain_old[_qp];
}
