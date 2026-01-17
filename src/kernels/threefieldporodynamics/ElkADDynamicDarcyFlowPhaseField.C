//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/**
 * ElkADDynamicDarcyFlowPhaseField - Dynamic Darcy flow kernel for phase-field damage coupling
 *
 * Implements: rho^f a^s + rho^f tau_t / phi * a^f + mu_f / kappa * v^f
 *
 * Newmark time integration:
 *   a^s_{n+1} = (u_{n+1} - u_n - dt*v^s_n) / (beta*dt^2) - a^s_n*(1-2*beta)/(2*beta)
 *   a^f_{n+1} = (w_{n+1} - w_n - dt*v^f_n) / (beta*dt^2) - a^f_n*(1-2*beta)/(2*beta)
 *   v^f_{n+1} = gamma*(w_{n+1} - w_n)/(beta*dt) + v^f_n*(1-gamma/beta) + a^f_n*dt*(1-gamma/(2*beta))
 */
#include "ElkADDynamicDarcyFlowPhaseField.h"

registerMooseObject("farmsApp", ElkADDynamicDarcyFlowPhaseField);

InputParameters
ElkADDynamicDarcyFlowPhaseField::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Dynamic Darcy flow kernel for three-field poroelastodynamics "
                             "with phase-field damage coupling. Uses RankTwoTensor permeability.");
  params.set<bool>("use_displaced_mesh") = false;
  params.addRequiredCoupledVar("skeletondisp", "Skeleton displacement variable");
  params.addRequiredCoupledVar("skeletonvel", "Skeleton velocity variable");
  params.addRequiredCoupledVar("skeletonaccel", "Skeleton acceleration variable");
  params.addRequiredCoupledVar("fluidvel", "Fluid relative velocity variable");
  params.addRequiredCoupledVar("fluidaccel", "Fluid relative acceleration variable");
  params.addRequiredParam<Real>("beta", "Newmark beta parameter");
  params.addRequiredParam<Real>("gamma", "Newmark gamma parameter");
  params.addRequiredParam<int>("component", "Component of displacement (0=x, 1=y, 2=z)");
  return params;
}

ElkADDynamicDarcyFlowPhaseField::ElkADDynamicDarcyFlowPhaseField(const InputParameters & parameters)
  : ADKernel(parameters),
    _rhof(getADMaterialProperty<Real>("rhof")),
    _nf(getADMaterialProperty<Real>("porosity")),
    _taut(getADMaterialProperty<Real>("tortosity")),
    _muf(getADMaterialProperty<Real>("viscosity")),
    _kappa(getADMaterialProperty<RankTwoTensor>("permeability")),
    _us(adCoupledValue("skeletondisp")),
    _us_old(coupledValueOld("skeletondisp")),
    _vs_old(coupledValueOld("skeletonvel")),
    _as_old(coupledValueOld("skeletonaccel")),
    _fluiddisp_old(valueOld()),
    _fluidvel_old(coupledValueOld("fluidvel")),
    _fluidaccel_old(coupledValueOld("fluidaccel")),
    _beta(getParam<Real>("beta")),
    _gamma(getParam<Real>("gamma")),
    _component(getParam<int>("component"))
{
}

ADReal
ElkADDynamicDarcyFlowPhaseField::computeQpResidual()
{
  // Neglect inertia term at t = 0
  if (_dt == 0)
    return 0.0;

  // Compute solid acceleration using Newmark method
  // a^s_{n+1} = (u_{n+1} - u_n - dt*v^s_n) / (beta*dt^2) - a^s_n*(1-2*beta)/(2*beta)
  ADReal as = (_us[_qp] - _us_old[_qp] - _dt * _vs_old[_qp]) / (_beta * _dt * _dt) -
              _as_old[_qp] * (1.0 - 2.0 * _beta) / (2.0 * _beta);

  // Compute fluid acceleration using Newmark method
  // a^f_{n+1} = (w_{n+1} - w_n - dt*v^f_n) / (beta*dt^2) - a^f_n*(1-2*beta)/(2*beta)
  ADReal af = (_u[_qp] - _fluiddisp_old[_qp] - _dt * _fluidvel_old[_qp]) / (_beta * _dt * _dt) -
              _fluidaccel_old[_qp] * (1.0 - 2.0 * _beta) / (2.0 * _beta);

  // Compute fluid velocity using Newmark method
  // v^f_{n+1} = gamma*(w_{n+1} - w_n)/(beta*dt) + v^f_n*(1-gamma/beta) + a^f_n*dt*(1-gamma/(2*beta))
  ADReal vf = _gamma * (_u[_qp] - _fluiddisp_old[_qp]) / (_beta * _dt) +
              _fluidvel_old[_qp] * (1.0 - _gamma / _beta) +
              _fluidaccel_old[_qp] * _dt * (1.0 - _gamma / (2.0 * _beta));

  // Compute residual: rho^f * a^s + rho^f * tau_t / phi * a^f + mu_f / kappa * v^f
  // Permeability is diagonal, use component-wise access
  return _test[_i][_qp] *
         (_rhof[_qp] * as + _rhof[_qp] * _taut[_qp] / _nf[_qp] * af +
          _muf[_qp] / _kappa[_qp](_component, _component) * vf);
}
