// //* This file is part of the MOOSE framework
// //* https://www.mooseframework.org
// //*
// //* All rights reserved, see COPYRIGHT for full restrictions
// //* https://github.com/idaholab/moose/blob/master/COPYRIGHT
// //*
// //* Licensed under LGPL 2.1, please see LICENSE for details
// //* https://www.gnu.org/licenses/lgpl-2.1.html

// #include "SemiPermeableFaultConditionBC.h"
// #include "Material.h" 

// registerMooseObject("MooseApp", SemiPermeableFaultConditionBC);

// InputParameters
// SemiPermeableFaultConditionBC::validParams()
// {
//   InputParameters params = DirichletBCBase::validParams();
//   params.addClassDescription("Imposes the essential boundary condition $u=g$, where $g$ "
//                              "is a constant, controllable value.");
//   params.addParam<std::string>("base_name", "Material property base name");
//   params.addParam<std::string>("side", "main", "Specifies which side to apply the condition: 'main' or 'secondary'.");
//   return params;
// }

// SemiPermeableFaultConditionBC::SemiPermeableFaultConditionBC(const InputParameters & parameters)
//   : DirichletBCBase(parameters),
//     _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
//     _across_flux_main(getMaterialPropertyByName<Real>(_base_name + "across_flux_main")),
//     _across_flux_sec(getMaterialPropertyByName<Real>(_base_name + "across_flux_sec")),
//     _side(getParam<std::string>("side"))
// {
// }

// Real
// SemiPermeableFaultConditionBC::computeQpValue()
// {
//   if (_side == "main")
//   {
//     return _across_flux_main[_qp];
//   }
//   else
//   {
//     return _across_flux_sec[_qp];
//   }
// }
