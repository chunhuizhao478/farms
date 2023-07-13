//Add constant value to variable

#pragma once

#include "AuxKernel.h"

class CompTJumpRateState : public AuxKernel
{
    public:

    static InputParameters validParams();
    CompTJumpRateState(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    const VariableValue & _coupled_val;

    Real _sliprate_bd;

};