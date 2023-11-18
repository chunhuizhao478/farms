/*
AuxKernel to calculate total prinpal strain rate 
*/

#pragma once

#include "AuxKernel.h"

class PrincipalStrainCalc : public AuxKernel
{
    public:

    static InputParameters validParams();
    PrincipalStrainCalc(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    //current value
    const VariableValue & _eps_total_11;
    const VariableValue & _eps_total_12;
    const VariableValue & _eps_total_22;

    //old value
    const VariableValue & _eps_total_11_old;
    const VariableValue & _eps_total_12_old;
    const VariableValue & _eps_total_22_old;

};