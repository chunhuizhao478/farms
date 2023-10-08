/*
Check Alpha and B within the range
*/

#pragma once

#include "AuxKernel.h"

class CheckAlphaBDev : public AuxKernel
{
    public:

    static InputParameters validParams();
    CheckAlphaBDev(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    const VariableValue & _coupled_val;

};