//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/**
 * Created by Chunhui Zhao, Jul 2nd, 2024
 * Implement rho^f a^s + rho^f tau_t / phi * a^f + mu_f / kappa * v^f terms
 * rho^f: fluid density (_rhof[_qp])
 * a^s: solid acceleration - Newmark method, a^s_{n+1} = ( u_{n+1} - u_{n} - \Delta t v^s_{n} ) / ( \beta \Delta t ^ 2 ) - a^s_{n} ( 1 - 2 \beta ) / ( 2 \beta )
 * tau_t : tortosity (_taut[_qp])
 * phi: porosity (_nf[_qp])
 * a^f: fluid acceleration - Newmark method, a^f_{n+1} = ( w_{n+1} - w_{n} - \Delta t v^f_{n} ) / ( \beta \Delta t ^ 2 ) - a^f_{n} ( 1 - 2 \beta ) / ( 2 \beta )
 * mu_f: viscosity (_muf[_qp])
 * kappa: permeability (_kappa[_qp])
 * v^f : darcy velocity, use Newmark Method 
 * v^f_{n+1} = \gamma ( u_{n+1} - u_{n} ) / ( \beta \Delta t ) + v^f_{n} ( 1 - \gamma / \beta ) + a^f_{n} \Delta t ( 1 - \gamma/(2 \beta))
 */
#include "ElkADDynamicDarcyFlowSmearedCracking.h"
#include "SubProblem.h"

registerMooseObject("farmsApp", ElkADDynamicDarcyFlowSmearedCracking);

InputParameters ElkADDynamicDarcyFlowSmearedCracking::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.set<bool>("use_displaced_mesh") = false;
    params.addRequiredCoupledVar("skeletondisp","skeleton displacement variable");
    params.addRequiredCoupledVar("skeletonvel","skeleton velocity variable");
    params.addRequiredCoupledVar("skeletonaccel","skeleton acceleration variable");
    params.addRequiredCoupledVar("fluidvel","fluid relative velocity variable");
    params.addRequiredCoupledVar("fluidaccel","fluid relative acceleration variable");
    params.addRequiredParam<Real>("beta","beta parameter");
    params.addRequiredParam<Real>("gamma","gamma parameter");
    params.addRequiredParam<int>("component","component of displacements");
    return params;
}

ElkADDynamicDarcyFlowSmearedCracking::ElkADDynamicDarcyFlowSmearedCracking(const InputParameters & parameters):ADKernel(parameters),
    _rhof(getADMaterialProperty<Real>("rhof")),
    _nf(getADMaterialProperty<Real>("porosity")),
    _taut(getADMaterialProperty<Real>("tortosity")),
    _muf(getADMaterialProperty<Real>("viscosity")),
    _kappa(getADMaterialProperty<RealTensorValue>("permeability")),
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
{}

ADReal
ElkADDynamicDarcyFlowSmearedCracking::computeQpResidual()
{
    //neglect inertia term at t = 0
    if (_dt == 0)
    return 0.0;

    //compute solid and fluid acceleration using Newmark method
    ADReal as = ( _us[_qp] -        _us_old[_qp] - _dt *       _vs_old[_qp] ) / ( _beta * _dt * _dt ) -         _as_old[_qp] * ( 1.0 - 2.0 * _beta ) / ( 2.0 * _beta );
    ADReal af = (  _u[_qp] - _fluiddisp_old[_qp] - _dt * _fluidvel_old[_qp] ) / ( _beta * _dt * _dt ) - _fluidaccel_old[_qp] * ( 1.0 - 2.0 * _beta ) / ( 2.0 * _beta );
    ADReal vf = _gamma*( _u[_qp] - _fluiddisp_old[_qp] ) / ( _beta * _dt ) + _fluidvel_old[_qp] * ( 1 - _gamma/_beta ) + _fluidaccel_old[_qp] * _dt * ( 1.0 - _gamma / ( 2.0 * _beta ));

    //compute terms
    //w is u (primary variable)
    //here for each direction, the _kappa is the global direction values (xx, yy, zz)
    return _test[_i][_qp]*( _rhof[_qp]*as + _rhof[_qp]*_taut[_qp]/_nf[_qp]*af + _muf[_qp]/_kappa[_qp](_component,_component)*vf ); 
}