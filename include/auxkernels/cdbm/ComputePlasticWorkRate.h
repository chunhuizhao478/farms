/*
AuxKernel of Computing Plastic Work Rate For each element
*/

#pragma once

#include "AuxKernel.h"

class ComputePlasticWorkRate : public AuxKernel
{
    public:

    static InputParameters validParams();
    ComputePlasticWorkRate(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    const VariableValue & _eps_p_11;
    const VariableValue & _eps_p_11_old;
    const VariableValue & _eps_p_12;
    const VariableValue & _eps_p_12_old;
    const VariableValue & _eps_p_22;
    const VariableValue & _eps_p_22_old;
    const VariableValue & _eps_p_33;
    const VariableValue & _eps_p_33_old;
    const VariableValue & _sts_change_11;
    const VariableValue & _sts_change_12;
    const VariableValue & _sts_change_22;
    const VariableValue & _sts_change_33;
    const VariableValue & _sts_initial_11;
    const VariableValue & _sts_initial_12;
    const VariableValue & _sts_initial_22;
    const VariableValue & _sts_initial_33;

};