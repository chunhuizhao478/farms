/*
Check Alpha and B within the range
*/

#pragma once

#include "AuxKernel.h"

class CheckAlphaB : public AuxKernel
{
    public:

    static InputParameters validParams();
    CheckAlphaB(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    const VariableValue & _coupled_val;

};