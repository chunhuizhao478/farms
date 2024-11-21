
#pragma once

#include "Kernel.h"
#include "PorousFlowDictator.h"

/**
 * Darcy advective flux for a fully-saturated,
 * single phase, single component fluid.
 * No upwinding or relative-permeability is used.
 */
class PorousFlowDarcyPressure : public Kernel
{
public:
  static InputParameters validParams();

  PorousFlowDarcyPressure(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  
  /// If true then the mobility contains the fluid density, otherwise it doesn't
  const bool _multiply_by_density;

  /// Fluid density for each phase (at the qp)
  const MaterialProperty<std::vector<Real>> & _density;

  /// Derivative of the fluid density for each phase wrt PorousFlow variables (at the qp)
  const MaterialProperty<std::vector<std::vector<Real>>> & _ddensity_dvar;

  /// Quadpoint pore pressure in each phase
  const MaterialProperty<std::vector<Real>> & _pp;

  // Derivative of porepressure in each phase wrt the PorousFlow variables
  const MaterialProperty<std::vector<std::vector<Real>>> & _dpp_dvar;

  /// PorousFlowDictator UserObject
  const PorousFlowDictator & _dictator;

  /// The spatial component
  const unsigned int _component;

};
