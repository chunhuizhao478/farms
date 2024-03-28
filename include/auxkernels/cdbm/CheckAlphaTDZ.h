/*
Check Alpha and B within the range
*/

#pragma once

#include "AuxKernel.h"

class CheckAlphaTDZ : public AuxKernel
{
    public:

    static InputParameters validParams();
    CheckAlphaTDZ(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    const VariableValue & _coupled_val;

    const VariableValue & _initial_damage;

};