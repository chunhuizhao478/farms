//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StateEvolutionKernel.h"
#include <cmath>

registerMooseObject("farmsApp", StateEvolutionKernel);

InputParameters
StateEvolutionKernel::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription(
      "Implements state evolution equation for rate-and-state friction. "
      "Aging law: dθ/dt = 1 - V*θ/Dc. "
      "Slip law: dθ/dt = -V*θ/Dc * ln(V*θ/Dc). "
      "Provides RHS contribution to be used with TimeDerivative kernel.");

  params.addRequiredCoupledVar("slip_rate", "The slip rate variable V (m/s)");
  params.addRequiredParam<Real>("Dc", "Critical slip distance Dc (m)");

  MooseEnum evolution_law("aging slip", "aging");
  params.addParam<MooseEnum>(
      "evolution_law", evolution_law, "State evolution law: 'aging' (default) or 'slip'");

  return params;
}

StateEvolutionKernel::StateEvolutionKernel(const InputParameters & parameters)
  : Kernel(parameters),
    _slip_rate(coupledValue("slip_rate")),
    _slip_rate_var(coupled("slip_rate")),
    _Dc(getParam<Real>("Dc")),
    _evolution_law(getParam<MooseEnum>("evolution_law"))
{
  if (_Dc <= 0.0)
    mooseError("Critical slip distance Dc must be positive");
}

Real
StateEvolutionKernel::computeQpResidual()
{
  // θ is the variable _u
  Real theta = _u[_qp];
  Real V = std::fabs(_slip_rate[_qp]);

  // Ensure positive state variable
  if (theta <= 0.0)
    theta = _Dc / std::max(V, 1e-20);

  Real rhs;
  if (_evolution_law == "aging")
  {
    // Aging law: dθ/dt = 1 - V*θ/Dc
    // Residual: -test * (1 - V*θ/Dc)
    rhs = 1.0 - V * theta / _Dc;
  }
  else // slip law
  {
    // Slip law: dθ/dt = -V*θ/Dc * ln(V*θ/Dc)
    Real Vtheta_Dc = V * theta / _Dc;
    if (Vtheta_Dc <= 0.0)
      Vtheta_Dc = 1e-20;
    rhs = -Vtheta_Dc * std::log(Vtheta_Dc);
  }

  // Negative sign because we're providing -RHS to balance TimeDerivative
  return -_test[_i][_qp] * rhs;
}

Real
StateEvolutionKernel::computeQpJacobian()
{
  // Derivative of residual w.r.t. θ (the variable)
  Real V = std::fabs(_slip_rate[_qp]);

  if (_evolution_law == "aging")
  {
    // d/dθ[-test * (1 - V*θ/Dc)] = -test * (-V/Dc) * φ = test * V/Dc * φ
    return _test[_i][_qp] * (V / _Dc) * _phi[_j][_qp];
  }
  else // slip law
  {
    // d/dθ[-test * (-V*θ/Dc * ln(V*θ/Dc))]
    // = test * V/Dc * (ln(V*θ/Dc) + 1) * φ
    Real theta = std::max(_u[_qp], _Dc / std::max(V, 1e-20));
    Real Vtheta_Dc = V * theta / _Dc;
    if (Vtheta_Dc <= 0.0)
      Vtheta_Dc = 1e-20;
    return _test[_i][_qp] * (V / _Dc) * (std::log(Vtheta_Dc) + 1.0) * _phi[_j][_qp];
  }
}

Real
StateEvolutionKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Off-diagonal Jacobian w.r.t. slip_rate
  // Since slip_rate is an AuxVariable (not solved), this is typically not needed
  // But included for completeness if using coupled variables
  if (jvar == _slip_rate_var)
  {
    Real theta = std::max(_u[_qp], 1e-20);

    if (_evolution_law == "aging")
    {
      // d/dV[-test * (1 - V*θ/Dc)] = test * θ/Dc * φ_V
      // But V is AuxVariable so φ_V = 0 in standard MOOSE
      return 0.0;
    }
    else // slip law
    {
      return 0.0;
    }
  }

  return 0.0;
}
