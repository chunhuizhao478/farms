/*
AuxKernel of Passing Variable
*/

#pragma once

#include "AuxKernel.h"

class CompXi3D : public AuxKernel
{
    public:

    static InputParameters validParams();
    CompXi3D(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    const VariableValue & _strainxx;
    const VariableValue & _strainxy;
    const VariableValue & _strainyy;
    const VariableValue & _strainxz;
    const VariableValue & _strainyz;
    const VariableValue & _strainzz;
};



