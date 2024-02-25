/*
Check Alpha and B within the range
*/

#pragma once

#include "AuxKernel.h"

class CheckAlphaDev : public AuxKernel
{
    public:

    static InputParameters validParams();
    CheckAlphaDev(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    Real _shear_modulus_o;
    Real _xi_0;
    Real _gamma_damaged_r;
    Real _factor;
    const VariableValue & _coupled_val;

};