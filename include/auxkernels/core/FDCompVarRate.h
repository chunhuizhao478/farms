//approximate rate using finite difference 

#pragma once

#include "AuxKernel.h"

class FDCompVarRate : public AuxKernel
{
    public:

    static InputParameters validParams();
    FDCompVarRate(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    const VariableValue & _coupled_val;
    const VariableValue & _coupled_val_old;

};