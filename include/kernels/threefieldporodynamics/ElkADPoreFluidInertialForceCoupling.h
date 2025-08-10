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
#pragma once

#include "ADKernel.h"
#include "Material.h"

//Forward Declarations

class ElkADPoreFluidInertialForceCoupling : public ADKernel
{
public:
    static InputParameters validParams();
    ElkADPoreFluidInertialForceCoupling(const InputParameters & parameters);
    
protected:
    virtual ADReal computeQpResidual() override;

private:
    const ADMaterialProperty<Real> & _rhof;  // fluid density
    const ADVariableValue & _fluiddisp;      // current fluid displacement
    const VariableValue   & _fluiddisp_old;  // old fluid displacement
    const VariableValue   & _fluidvel_old;   // old fluid velocity
    const VariableValue   & _fluidaccel_old; // old fluid acceleration
    ADReal _beta;
};
