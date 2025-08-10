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
 * Implement rho^f a^f term
 * rho^f: fluid density
 * a^f: fluid acceleration
 * Nemark method gives a^f_{n+1} = ( w_{n+1} - w_{n} - \Delta t v^f_{n} ) / ( \beta \Delta t ^ 2 ) - a^f_{n} ( 1 - 2 \beta ) / ( 2 \beta )
 */
#include "ElkADPoreFluidInertialForceCoupling.h"
#include "SubProblem.h"

registerMooseObject("farmsApp", ElkADPoreFluidInertialForceCoupling);

InputParameters ElkADPoreFluidInertialForceCoupling::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.set<bool>("use_displaced_mesh") = false;
    params.addRequiredCoupledVar("fluiddisp","fluid displacement variable");
    params.addRequiredCoupledVar("fluidvel","fluid velocity variable");
    params.addRequiredCoupledVar("fluidaccel","fluid relative acceleration variable");
    params.addRequiredParam<Real>("beta","beta parameter");
    return params;
}

ElkADPoreFluidInertialForceCoupling::ElkADPoreFluidInertialForceCoupling(const InputParameters &
parameters):ADKernel(parameters),
    _rhof(getADMaterialProperty<Real>("rhof")),
    _fluiddisp(adCoupledValue("fluiddisp")),
    _fluiddisp_old(coupledValueOld("fluiddisp")),
    _fluidvel_old(coupledValueOld("fluidvel")),
    _fluidaccel_old(coupledValueOld("fluidaccel")),
    _beta(getParam<Real>("beta"))
{}

ADReal
ElkADPoreFluidInertialForceCoupling::computeQpResidual()
{
    //neglect inertia term at t = 0
    if (_dt == 0)
    return 0.0;

    //compute fluid acceleration
    ADReal af = ( _fluiddisp[_qp] - _fluiddisp_old[_qp] - _dt * _fluidvel_old[_qp] ) / ( _beta * _dt * _dt ) - _fluidaccel_old[_qp] * ( 1.0 - 2.0 * _beta ) / ( 2.0 * _beta );
    
    return _test[_i][_qp]*_rhof[_qp]*af;
}