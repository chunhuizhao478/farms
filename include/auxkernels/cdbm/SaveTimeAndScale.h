/*
Check Alpha and B within the range
*/

#pragma once

#include "AuxKernel.h"

class SaveTimeAndScale : public AuxKernel
{
    public:

    static InputParameters validParams();
    SaveTimeAndScale(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

};