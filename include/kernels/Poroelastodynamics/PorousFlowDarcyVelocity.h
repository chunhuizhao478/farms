

#pragma once

#include "Kernel.h"
#include "PorousFlowDictator.h"
#include "RankTwoTensor.h"


/**
 * Darcy advective flux for a fully-saturated,
 */
class PorousFlowDarcyVelocity : public Kernel
{
public:
  static InputParameters validParams();

  PorousFlowDarcyVelocity(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /**
   * The mobility of the fluid = density / viscosity
   */
  virtual Real mobility() const;

  /**
   * The derivative of the mobility with respect to the PorousFlow variable pvar
   * @param pvar Take the derivative with respect to this PorousFlow variable
   */
  virtual Real dmobility(unsigned pvar) const;

  /// If true then the mobility contains the fluid density, otherwise it doesn't
  const bool _multiply_by_density;

  /// Permeability of porous material
  const MaterialProperty<RealTensorValue> & _permeability;

  /// d(permeabiity)/d(PorousFlow variable)
  const MaterialProperty<std::vector<RealTensorValue>> & _dpermeability_dvar;

  /// d(permeabiity)/d(grad(PorousFlow variable))
  const MaterialProperty<std::vector<std::vector<RealTensorValue>>> & _dpermeability_dgradvar;

  /// Fluid density for each phase (at the qp)
  const MaterialProperty<std::vector<Real>> & _density;

  /// Derivative of the fluid density for each phase wrt PorousFlow variables (at the qp)
  const MaterialProperty<std::vector<std::vector<Real>>> & _ddensity_dvar;

  /// Viscosity of the fluid at the qp
  const MaterialProperty<std::vector<Real>> & _viscosity;

  /// Derivative of the fluid viscosity  wrt PorousFlow variables
  const MaterialProperty<std::vector<std::vector<Real>>> & _dviscosity_dvar;

  /// PorousFlowDictator UserObject
  const PorousFlowDictator & _dictator;
  
   /// An integer corresponding to the direction this kernel acts in
  const unsigned int _component;

  /// Flag to check whether permeabiity derivatives are non-zero
  const bool _perm_derivs;
};