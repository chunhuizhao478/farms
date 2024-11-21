// //* This file is part of the MOOSE framework
// //* https://www.mooseframework.org
// //*
// //* All rights reserved, see COPYRIGHT for full restrictions
// //* https://github.com/idaholab/moose/blob/master/COPYRIGHT
// //*
// //* Licensed under LGPL 2.1, please see LICENSE for details
// //* https://www.gnu.org/licenses/lgpl-2.1.html

// #pragma once

// #include "DirichletBCBase.h"

// /**
//  * Boundary condition of a Dirichlet type
//  *
//  * Sets the value in the node
//  */
// class SemiPermeableFaultConditionBC : public DirichletBCBase
// {
// public:
//   static InputParameters validParams();

//   SemiPermeableFaultConditionBC(const InputParameters & parameters);

// protected:
//   virtual Real computeQpValue() override;

//   const std::string _base_name;
//   const MaterialProperty<Real> & _across_flux_main;
//   const MaterialProperty<Real> & _across_flux_sec;
//   std::string _side;  // Declare the _side member variable
// };