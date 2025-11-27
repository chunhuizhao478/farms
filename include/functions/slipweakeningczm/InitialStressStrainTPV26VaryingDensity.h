// This file declares the TPV26 initial stress/strain Function with depth-varying density (TPV32 profile).
#pragma once

#include "Function.h"

class InitialStressStrainTPV26VaryingDensity : public Function
{
public:
	static InputParameters validParams();
	InitialStressStrainTPV26VaryingDensity(const InputParameters & parameters);

	Real value(Real t, const Point & p) const override;

private:
	// TPV32 piecewise-linear density (kg/m^3) at depth d (m)
	Real tpv32_density_at_depth(Real d) const;

	// Inputs/flags
	const Real _i;
	const Real _j;
	const Real _lambda_o;
	const Real _shear_modulus_o;
	const Real _fluid_density;
	const Real _gravity;
	const Real _bxx;
	const Real _byy;
	const Real _bxy;

	const bool _get_initial_stress;
	const bool _get_initial_strain;
	const bool _get_shear_overstress;
	const bool _get_fluid_pressure;

	const bool _use_tapering;
	const Real _tapering_depth_A;
	const Real _tapering_depth_B;

	const bool _use_overpressure;
	const Real _overpressure_depth_A;
	const Real _overpressure_depth_B;
	const Real _overpressure_rho_ref;

	const bool _flip_sign;
	const Real _depth_offset;
};
