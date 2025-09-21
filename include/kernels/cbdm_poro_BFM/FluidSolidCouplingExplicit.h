//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

/// Fluid-solid coupling kernel for poromechanical problems
///
/// This kernel solves the pore pressure equation with coupling to displacement fields.
/// It computes the volumetric strain contribution from elastic deformation and couples
/// it to the fluid flow equation.
///
class FluidSolidCouplingExplicit : public Kernel
{
public:
  static InputParameters validParams();
  FluidSolidCouplingExplicit(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Number of displacement components
  const unsigned int _ndisp;

  /// Variable numbers for displacement variables (for off-diagonal Jacobian)
  std::vector<unsigned int> _disp_var;

  /// First elastic strain invariant from material (current)
  const MaterialProperty<Real> & _I1;
  
  /// Old first elastic strain invariant from material  
  const MaterialProperty<Real> & _I1_older;

/// Deformation gradient and plastic deformation gradient (needed for Jacobian)
  const MaterialProperty<RankTwoTensor> & _F;
  const MaterialProperty<RankTwoTensor> & _dI1dF;

  /// Poroelastic properties
  const MaterialProperty<Real> & _fluid_solid_coupling;
};