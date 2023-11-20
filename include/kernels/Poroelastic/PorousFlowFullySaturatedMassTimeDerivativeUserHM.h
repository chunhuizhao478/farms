#pragma once

#include "TimeKernel.h"
#include "PorousFlowDictator.h"

/**
 * Time derivative of fluid mass suitable for fully-saturated,
 * single-phase, single-component simulations.
 * No mass-lumping is used
 */
class PorousFlowFullySaturatedMassTimeDerivativeUserHM : public TimeKernel
{
public:
  static InputParameters validParams();

  PorousFlowFullySaturatedMassTimeDerivativeUserHM (const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// Jacobian contribution for the PorousFlow variable pvar
  Real computeQpJac(unsigned int pvar);

  /// PorousFlowDictator UserObject
  const PorousFlowDictator & _dictator;

  /// Whether the Variable for this Kernel is a PorousFlow variable
  const bool _var_is_porflow_var;

  /// If true then the Kernel is the time derivative of the fluid mass, otherwise it is the derivative of the fluid volume
  const bool _multiply_by_density;

  /// Determines whether mechanical and/or thermal contributions should be added to the residual
  const enum class CouplingTypeEnum {
    HydroMechanical
  } _coupling_type;

  /// Whether mechanical contributions should be added to the residual
  const bool _includes_mechanical;

  /// Biot coefficient (used in simulations involving Mechanical deformations)
  const Real _biot_coefficient;

  /// Constant Biot modulus
  const MaterialProperty<Real> & _biot_modulus;

  /// Quadpoint fluid density for each phase
  const MaterialProperty<std::vector<Real>> * const _fluid_density;

  /// derivative of fluid density for each phase with respect to the PorousFlow variables
  const MaterialProperty<std::vector<std::vector<Real>>> * const _dfluid_density_dvar;

  /// Quadpoint pore pressure in each phase
  const MaterialProperty<std::vector<Real>> & _pp;

  /// Old value of quadpoint pore pressure in each phase
  const MaterialProperty<std::vector<Real>> & _pp_old;

  /// Derivative of porepressure in each phase wrt the PorousFlow variables
  const MaterialProperty<std::vector<std::vector<Real>>> & _dpp_dvar;

};