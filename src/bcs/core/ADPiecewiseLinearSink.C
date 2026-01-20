//* This file is part of the FARMS application
//*
//* AD version of PorousFlowPiecewiseLinearSink from MOOSE PorousFlow module
//* Simplified for three-field poroelastodynamics without PorousFlow dependency

#include "ADPiecewiseLinearSink.h"
#include "Function.h"

registerMooseObject("farmsApp", ADPiecewiseLinearSink);

InputParameters
ADPiecewiseLinearSink::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription(
      "Applies a flux sink to a boundary. The flux is computed as: "
      "flux = C * g(p - PT_shift), where C is flux_function (conductance), "
      "g is the piecewise linear function defined by pt_vals and multipliers, "
      "and PT_shift is the reference pressure. "
      "Positive flux is OUT of the domain (sink), negative is IN (source). "
      "This is an AD version for use with three-field poroelastodynamics.");

  params.addRequiredParam<std::vector<Real>>(
      "pt_vals",
      "Tuple of x-coordinates defining the piecewise linear function g(x). "
      "Must be monotonically increasing. For a simple Robin BC where g(x)=x, "
      "use pt_vals='-1e9 1e9' with multipliers='-1e9 1e9'.");
  params.addRequiredParam<std::vector<Real>>(
      "multipliers",
      "Tuple of y-coordinates defining the piecewise linear function g(x). "
      "Must have the same size as pt_vals.");

  params.addParam<FunctionName>(
      "flux_function",
      1.0,
      "The conductance C (can be a Function of time and space, or a constant). "
      "The flux is: flux = C * g(p - PT_shift). "
      "Units should be consistent with your problem (e.g., m/s/Pa for pressure BC).");

  params.addCoupledVar("pressure_variable",
                       "The pore pressure variable. If not supplied, the BC variable is used.");

  params.addParam<FunctionName>("PT_shift",
                                "0",
                                "Reference/boundary pressure function. Can be a constant or "
                                "a Function of time and space. The piecewise linear function "
                                "is evaluated at (p - PT_shift). For a Robin BC with reference "
                                "pressure p_ref(t,x), set PT_shift to that function.");

  params.addParam<unsigned int>("fluid_phase",
                                0,
                                "The fluid phase number (for compatibility with PorousFlow). "
                                "Not used in this simplified implementation.");

  params.addParam<bool>("use_mobility",
                        false,
                        "If true, flux is multiplied by (permeability_nn * density / viscosity). "
                        "Requires permeability, fluid_density, and fluid_viscosity materials.");

  params.addParam<MaterialPropertyName>("permeability",
                                        "permeability",
                                        "The permeability tensor material property name.");
  params.addParam<MaterialPropertyName>("fluid_density",
                                        "fluid_density",
                                        "The fluid density material property name.");
  params.addParam<MaterialPropertyName>("fluid_viscosity",
                                        "fluid_viscosity",
                                        "The fluid viscosity material property name.");

  return params;
}

ADPiecewiseLinearSink::ADPiecewiseLinearSink(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _flux_function(getFunction("flux_function")),
    _sink_func(getParam<std::vector<Real>>("pt_vals"),
               getParam<std::vector<Real>>("multipliers")),
    _pressure(isParamValid("pressure_variable") ? adCoupledValue("pressure_variable") : _u),
    _PT_shift_function(getFunction("PT_shift")),
    _fluid_phase(getParam<unsigned int>("fluid_phase")),
    _use_mobility(getParam<bool>("use_mobility")),
    _permeability(_use_mobility ? &getADMaterialProperty<RankTwoTensor>("permeability") : nullptr),
    _fluid_density(_use_mobility ? &getADMaterialProperty<Real>("fluid_density") : nullptr),
    _fluid_viscosity(_use_mobility ? &getADMaterialProperty<Real>("fluid_viscosity") : nullptr)
{
  // Validate pt_vals and multipliers have the same size
  const auto & pt_vals = getParam<std::vector<Real>>("pt_vals");
  const auto & multipliers = getParam<std::vector<Real>>("multipliers");

  if (pt_vals.size() != multipliers.size())
    paramError("multipliers",
               "The 'pt_vals' and 'multipliers' vectors must have the same size. "
               "pt_vals has ",
               pt_vals.size(),
               " entries and multipliers has ",
               multipliers.size(),
               " entries.");

  if (pt_vals.size() < 2)
    paramError("pt_vals",
               "The 'pt_vals' vector must have at least 2 entries for interpolation.");

  // Check that pt_vals is monotonically increasing
  for (unsigned int i = 1; i < pt_vals.size(); ++i)
    if (pt_vals[i] <= pt_vals[i - 1])
      paramError("pt_vals",
                 "The 'pt_vals' vector must be monotonically increasing. "
                 "Entry ",
                 i,
                 " (",
                 pt_vals[i],
                 ") is not greater than entry ",
                 i - 1,
                 " (",
                 pt_vals[i - 1],
                 ").");

  // Error checking for mobility
  if (_use_mobility)
  {
    if (!hasMaterialProperty<RankTwoTensor>("permeability"))
      paramError("use_mobility",
                 "use_mobility=true requires a 'permeability' material property.");
    if (!hasMaterialProperty<Real>("fluid_density"))
      paramError("use_mobility",
                 "use_mobility=true requires a 'fluid_density' material property.");
    if (!hasMaterialProperty<Real>("fluid_viscosity"))
      paramError("use_mobility",
                 "use_mobility=true requires a 'fluid_viscosity' material property.");
  }
}

ADReal
ADPiecewiseLinearSink::computeQpResidual()
{
  // Get reference pressure from function (can be time/space varying)
  const Real PT_shift = _PT_shift_function.value(_t, _q_point[_qp]);

  // Get pressure value and apply shift: x = p - PT_shift
  const ADReal x = _pressure[_qp] - PT_shift;

  // Evaluate the piecewise linear function g(x) at the shifted pressure
  // Note: LinearInterpolation::sample returns Real, not ADReal
  // We need to handle the AD derivative manually using first-order Taylor expansion
  const Real x_val = MetaPhysicL::raw_value(x);
  const Real g_val = _sink_func.sample(x_val);
  const Real dg_dx = _sink_func.sampleDerivative(x_val);

  // Compute the flux: flux = C * g(p - PT_shift)
  // The conductance C from flux_function
  const Real C = _flux_function.value(_t, _q_point[_qp]);

  // Flux with AD sensitivity through pressure
  // Using Taylor expansion: g(x) ≈ g(x_val) + dg/dx * (x - x_val)
  ADReal flux = C * (g_val + dg_dx * (x - x_val));

  // Apply mobility if requested: flux *= k_nn * rho / mu
  if (_use_mobility)
  {
    // Project permeability to normal direction: k_nn = n^T * K * n
    const ADReal k_nn = (_normals[_qp] * ((*_permeability)[_qp] * _normals[_qp]));
    const ADReal mobility = (*_fluid_density)[_qp] * k_nn / (*_fluid_viscosity)[_qp];
    flux *= mobility;
  }

  // Return the residual contribution: test * flux
  // Positive flux means mass leaving the domain (sink)
  return _test[_i][_qp] * flux;
}
