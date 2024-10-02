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

    const MaterialProperty<RankTwoTensor> & _mechanical_strain;
};



