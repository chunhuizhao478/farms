//approximate rate using finite difference 

#pragma once

#include "AuxKernel.h"

class FDCompVarRate2 : public AuxKernel
{
    public:

    static InputParameters validParams();
    FDCompVarRate2(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    const VariableValue & _coupled_val;
    const VariableValue & _coupled_val_old;

};