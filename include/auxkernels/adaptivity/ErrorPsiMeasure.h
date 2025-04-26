/*
Compute error measure of strain energy density for phase-field between two AM iterations
*/

#pragma once

#include "AuxKernel.h"

class ErrorPsiMeasure : public AuxKernel
{
    public:

    static InputParameters validParams();
    ErrorPsiMeasure(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    const ADMaterialProperty<Real> & _psie_active;
    const MaterialProperty<Real> & _psie_active_old;

};