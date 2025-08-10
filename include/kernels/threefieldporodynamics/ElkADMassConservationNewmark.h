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
 * Implement div v^s term
 * v^s: velocity of solid skeleton
 * Newmark method, see section 2.4.2.1 equation (52) (53) 
 */
#pragma once

#include "ADKernel.h"

class ElkADMassConservationNewmark : public ADKernel
{
public:
    static InputParameters validParams();
    ElkADMassConservationNewmark(const InputParameters & parameters);

protected:
    virtual ADReal computeQpResidual() override;

private:
    const ADVariableGradient & _grad_ux; // gradient of displacement
    const ADVariableGradient & _grad_uy;
    const ADVariableGradient & _grad_uz;
    const VariableGradient & _grad_ux_old; // gradient of previous displacement
    const VariableGradient & _grad_uy_old;
    const VariableGradient & _grad_uz_old;
    const VariableGradient & _grad_vx_old; // gradient of previous velocity
    const VariableGradient & _grad_vy_old;
    const VariableGradient & _grad_vz_old;
    const VariableGradient & _grad_ax_old; // gradient of previous acceleration
    const VariableGradient & _grad_ay_old;
    const VariableGradient & _grad_az_old;
    const ADMaterialProperty<Real> & _biot_alpha; //biot coefficient
    ADReal _beta;
    ADReal _gamma;
    /// whether or not multiply biot coefficient by weak form
    const bool _multiply_biot_coefficient;
    const ADMaterialProperty<Real> & _biot_modulus; //biot modulus
};