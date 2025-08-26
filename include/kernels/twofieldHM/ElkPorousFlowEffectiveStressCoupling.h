//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"
#include "PorousFlowDictator.h"

/**
 * ElkPorousFlowEffectiveStressCoupling computes
 * -alpha*effective_porepressure*grad_component(test)
 * where alpha is either a constant parameter or a material property
 * named "biot_coefficient" when use_damaged_biot=true.
 */
class ElkPorousFlowEffectiveStressCoupling : public Kernel
{
public:
  static InputParameters validParams();

  ElkPorousFlowEffectiveStressCoupling(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// The PorousFlow dictator that holds global info about the simulation
  const PorousFlowDictator & _dictator;

  /// Constant Biot coefficient (used if !_use_damaged_biot)
  const Real _coefficient_const;

  /// Whether to read Biot coefficient from material property
  const bool _use_damaged_biot;

  /// Material property Biot coefficient (valid if _use_damaged_biot)
  const MaterialProperty<Real> * _biot_coeff_mp;

  /// The spatial component
  const unsigned int _component;

  /// Effective porepressure
  const MaterialProperty<Real> & _pf;

  /// d(effective porepressure)/(d porflow variable)
  const MaterialProperty<std::vector<Real>> & _dpf_dvar;

  /// Whether an RZ coordinate system is being used
  const bool _rz;

  /// Helper to get alpha at this qp
  inline Real biot() const { return _use_damaged_biot ? (*_biot_coeff_mp)[_qp] : _coefficient_const; }
};
