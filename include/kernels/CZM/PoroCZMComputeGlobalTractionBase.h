//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceMaterial.h"
/**
 * Base class traction computing the traction used to impose equilibrium and its derivatives  w.r.t.
 * the global displacement jump, starting from the values provided from any CZM constituive material
 * model.
 */
class PoroCZMComputeGlobalTractionBase : public InterfaceMaterial
{
public:
  static InputParameters validParams();
  PoroCZMComputeGlobalTractionBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// method computing the equilibrium traction and its derivatives
  virtual void computeEquilibriumTractionAndPressureAndDerivatives() = 0;

  /// Base name of the material system
  const std::string _base_name;

  /// the value of the traction in global and interface coordinates
  ///@{
  MaterialProperty<RealVectorValue> & _traction_global;
  const MaterialProperty<RealVectorValue> & _interface_traction;
  ///@}

  /// the traction's derivatives w.r.t. the displacement and fluid velocity jump in global and interface coordinates
  ///@{
  MaterialProperty<RankTwoTensor> & _dtraction_djump_global;
  MaterialProperty<RankTwoTensor> & _dtraction_djump_global_vf;
  MaterialProperty<RealVectorValue> & _dtraction_dpressure_global;
  const MaterialProperty<RankTwoTensor> & _dinterface_traction_djump;
  const MaterialProperty<RankTwoTensor> & _dinterface_traction_djump_vf;
  const MaterialProperty<RealVectorValue> & _dinterface_traction_dpressure;
  ///@}

  /// the value of the traction in global and interface coordinates
  ///@{
  MaterialProperty<Real> & _pressure_global;
  const MaterialProperty<Real> & _interface_pressure;
  ///@}

  /// the traction's derivatives w.r.t. the displacement jump in global and interface coordinates
  ///@{
  MaterialProperty<RealVectorValue> & _dpressure_djump_global;
  MaterialProperty<RealVectorValue> & _dpressure_djump_global_vf;
  const MaterialProperty<RealVectorValue> & _dinterface_pressure_djump;
  const MaterialProperty<RealVectorValue> & _dinterface_pressure_djump_vf;
  ///@}

  /// the rotation matrix trnasforming from interface to global coordinates
  const MaterialProperty<RankTwoTensor> & _czm_total_rotation;
};