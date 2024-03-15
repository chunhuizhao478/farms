/*
AuxKernel of Passing Variable
*/

#pragma once

#include "AuxKernel.h"

class CompXi : public AuxKernel
{
    public:

    static InputParameters validParams();
    CompXi(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    const VariableValue & _strainxx;
    const VariableValue & _strainxy;
    const VariableValue & _strainyy;
};



