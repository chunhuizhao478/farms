//* This file is part of the FARMS application
//*
//* AD version of PorousFlowPiecewiseLinearSink from MOOSE PorousFlow module
//* Simplified for three-field poroelastodynamics without PorousFlow dependency

#pragma once

#include "ADIntegratedBC.h"
#include "LinearInterpolation.h"

/**
 * ADPiecewiseLinearSink applies a flux sink to a boundary.
 * The flux is computed as:
 *   flux = C * g(p - PT_shift)
 *
 * where:
 *   - C is the flux_function (conductance, can be a Function of time and space)
 *   - g is the piecewise linear function defined by pt_vals and multipliers
 *   - p is the pore pressure (either the variable or a coupled variable)
 *   - PT_shift is the reference/boundary pressure
 *
 * Example: For a simple Robin BC with conductance C and reference pressure p_ref:
 *   pt_vals = '-1e9 1e9'
 *   multipliers = '-1e9 1e9'  (this gives g(x) = x)
 *   PT_shift = p_ref
 *   flux_function = C
 *   => flux = C * (p - p_ref)
 *
 * Positive flux is OUT of the domain (sink), negative is IN (source).
 *
 * This BC is useful for:
 *   - Robin-type pressure boundaries
 *   - Semi-permeable boundaries where flux depends on pressure difference
 *   - Pressure-dependent leakage
 */
class ADPiecewiseLinearSink : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ADPiecewiseLinearSink(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// The flux function / conductance (can vary in time and space, or be constant)
  const Function & _flux_function;

  /// Piecewise-linear function g(x) that multiplies the flux based on pressure
  const LinearInterpolation _sink_func;

  /// The pore pressure variable value at quadrature points
  const ADVariableValue & _pressure;

  /// Reference pressure function (can be constant or time/space varying)
  const Function & _PT_shift_function;

  /// Fluid phase number (for compatibility, not used in simplified version)
  const unsigned int _fluid_phase;

  /// Whether to use mobility (permeability * density / viscosity)
  const bool _use_mobility;

  /// Permeability (tensor) - optional
  const ADMaterialProperty<RankTwoTensor> * const _permeability;

  /// Fluid density - optional
  const ADMaterialProperty<Real> * const _fluid_density;

  /// Fluid viscosity - optional
  const ADMaterialProperty<Real> * const _fluid_viscosity;
};
