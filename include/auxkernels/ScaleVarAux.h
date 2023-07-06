#pragma once

#include "AuxKernel.h"

class ScaleVarAux : public AuxKernel
{
    public:

    static InputParameters validParams();
    ScaleVarAux(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    const VariableValue & _coupled_val;

    const VariableValue & _scale;

};