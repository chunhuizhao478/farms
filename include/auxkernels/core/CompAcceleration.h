/*
AuxKernel of Passing Variable Time Derivative 
*/

#pragma once

#include "AuxKernel.h"

class CompAcceleration : public AuxKernel
{
    public:

    static InputParameters validParams();
    CompAcceleration(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    const VariableValue & _coupled_val;

};



