/*
AuxKernel to track the Cd evolution based on strain rate
*/

#pragma once

#include "AuxKernel.h"

class TrackCdEvoAux : public AuxKernel
{
    public:

    static InputParameters validParams();
    TrackCdEvoAux(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    Real _Cd_min; //minimum Cd value for small strain (e < 1e-4)
    /// threshold mechanical strain (Cd remains constant below this value)
    Real _mechanical_strain_rate_threshold;
    /// mechanical strain rate
    const VariableValue & _mechanical_strain_rate;
    /// power index
    Real _m;
    Real _scale;
    int _option;
    Real _Cd_constant;

};