//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SmallStrainFluidSolidCouplingExplicit.h"

registerMooseObject("farmsApp", SmallStrainFluidSolidCouplingExplicit);

InputParameters
SmallStrainFluidSolidCouplingExplicit::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Fluid-solid coupling kernel for poromechanical problems");
  
  params.addRequiredCoupledVar("displacements", "The displacement variables");
  params.addParam<std::string>("base_name", "", "Optional parameter for multiple material systems");
  
  return params;
}

SmallStrainFluidSolidCouplingExplicit::SmallStrainFluidSolidCouplingExplicit(const InputParameters & parameters)
  : Kernel(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _I1(getMaterialProperty<Real>(getParam<std::string>("base_name") + "I1")),
    _I1_older(getMaterialPropertyOlder<Real>(getParam<std::string>("base_name") + "I1")),
    _fluid_solid_coupling(getMaterialProperty<Real>("fluid_solid_coupling"))
{
  // Get displacement variable numbers for off-diagonal Jacobian
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
SmallStrainFluidSolidCouplingExplicit::computeQpResidual()
{
  // First Invariant (I1): The trace of the strain tensor
  Real I1 = _I1[_qp];           // Current time step
  Real I1_older = _I1_older[_qp];   // Previous time step

  // Compute the volumetric coupling contribution
  // This represents: (term11/term33) * (I1 - I1_old) / dt
  Real volumetric_rate = _fluid_solid_coupling[_qp] * (I1 - I1_older) / _dt / 2;
  
  // The residual contribution 
  // R = ∫ φ * volumetric_rate dV
  return _test[_i][_qp] * volumetric_rate;
}

Real
SmallStrainFluidSolidCouplingExplicit::computeQpJacobian()
{
  return 0.0;
}

Real 
SmallStrainFluidSolidCouplingExplicit::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    if (jvar == _disp_var[i])
    {
      // ∂R/∂I1 = test * Kf / dt
      Real dR_dI1 = _test[_i][_qp] * _fluid_solid_coupling[_qp] / _dt / 2;

      // For small strain: I1 = ∇·u = ∂u_x/∂x + ∂u_y/∂y + ∂u_z/∂z
      // dI1/du_i = ∂(∂u_i/∂x_i)/∂u_i = ∂φ_j/∂x_i
      // The derivative with respect to displacement component i 
      // is just the i-th component of the shape function gradient
      Real dI1_du = _grad_phi[_j][_qp](i);
      
      return dR_dI1 * dI1_du;
    }
  }

  return 0.0;
}
