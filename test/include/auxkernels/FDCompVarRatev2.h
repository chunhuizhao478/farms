//approximate rate using finite difference 

#pragma once

#include "AuxKernel.h"

class FDCompVarRatev2 : public AuxKernel
{
    public:

    static InputParameters validParams();
    FDCompVarRatev2(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    const VariableValue & _coupled_val;
    const VariableValue & _coupled_val_old;

};