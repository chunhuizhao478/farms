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
#pragma once

#include "ADKernel.h"
#include "Material.h"

//Forward Declarations

class ElkADDynamicDarcyFlowSmearedCracking : public ADKernel
{
public:
    static InputParameters validParams();
    ElkADDynamicDarcyFlowSmearedCracking(const InputParameters & parameters);

protected:
    virtual ADReal computeQpResidual() override;

private:
    const ADMaterialProperty<Real> & _rhof;  // fluid density
    const ADMaterialProperty<Real> & _nf;    // porosity
    const ADMaterialProperty<Real> & _taut;  // tortosity
    const ADMaterialProperty<Real> & _muf;   // viscosity
    const ADMaterialProperty<RealTensorValue> & _kappa; // permeability //This is tensor

    const ADVariableValue & _us;             // current skeleton displacement
    const VariableValue   & _us_old;         // old skeleton displacement
    const VariableValue   & _vs_old;         // old skeleton velocity
    const VariableValue   & _as_old;         // old skeleton acceleration
    
    const VariableValue   & _fluiddisp_old;  // old fluid displacement
    const VariableValue   & _fluidvel_old;   // old fluid velocity
    const VariableValue   & _fluidaccel_old; // old fluid acceleration
    
    ADReal _beta;
    ADReal _gamma;

    /// choose an equation
    const int _component;
};
