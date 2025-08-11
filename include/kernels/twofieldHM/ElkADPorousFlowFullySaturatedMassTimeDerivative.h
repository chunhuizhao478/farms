//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/**
 * Created by Chunhui Zhao, Aug 11th, 2025
 * Implement alpha * volumetric strain rate * test + 1 / biot modulus * pressure rate * test
 * alpha: biot coefficient
 */
#pragma once

#include "ADTimeKernelValue.h"
#include "Material.h"

//Forward Declarations

class ElkADPorousFlowFullySaturatedMassTimeDerivative : public ADTimeKernelValue
{
public:
    static InputParameters validParams();
    ElkADPorousFlowFullySaturatedMassTimeDerivative(const InputParameters & parameters);
    
protected:
    // Return the scalar that multiplies the test function at this qp
    virtual ADReal precomputeQpResidual() override;

private:
    const ADMaterialProperty<Real> & _biot_modulus; //biot modulus
    const ADMaterialProperty<Real> & _biot_coefficient; //biot coefficient
    const ADMaterialProperty<Real> & _vol_strain; //volumetric strain
    const MaterialProperty<Real> & _vol_strain_old; //old volumetric strain
};