#pragma once

#include "AuxKernel.h"

class CompDispRateState : public AuxKernel
{
    public:

    static InputParameters validParams();
    CompDispRateState(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    //current displacement change measured from initial
    const VariableValue & _current_disp;
    //background velocity (add each time step)
    const VariableValue & _vel_const;

};