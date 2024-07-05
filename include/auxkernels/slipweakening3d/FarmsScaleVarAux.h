#pragma once

#include "AuxKernel.h"

class FarmsScaleVarAux : public AuxKernel
{
    public:

    static InputParameters validParams();
    FarmsScaleVarAux(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    const VariableValue & _coupled_val;

    /// The volume (or length) of the current element
    const Real & _current_elem_volume;

    /// The volume (or length) of the current side
    const Real & _current_side_volume;

};