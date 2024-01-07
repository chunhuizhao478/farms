/*
AuxKernel of Passing Variable
*/

#pragma once

#include "AuxKernel.h"

class CompArctanyx : public AuxKernel
{
    public:

    static InputParameters validParams();
    CompArctanyx(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

};



