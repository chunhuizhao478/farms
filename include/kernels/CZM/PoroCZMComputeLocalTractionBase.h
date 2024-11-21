
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
 * Base class used to implement traction separetion laws. All traction separetion laws
 * shall compute the interface traction using the interface coordinate system, and
 * effective traction derivtives w.r.t. to the interface displacement jump and interface fluid veclocity jump. 
 * Interface traction and related derivatives should be implemented overriding the computeInterfaceTractionAndDerivatives method.
 * The interface coordinate system assumes the three component of the traction and
 * disaplcement and fluid velocity jump being ordered as [N,S1,S2], where N is the normal component and S1, S2 two
 * orthogonal tangential components. The model also assumes isotropic behavior in the tangential
 * directions.
 */
class PoroCZMComputeLocalTractionBase : public InterfaceMaterial
{
public:
  static InputParameters validParams();
  PoroCZMComputeLocalTractionBase(const InputParameters & parameters);

protected:
  void initQpStatefulProperties() override;
  void computeQpProperties() override;

  /// Compute the local traction and derivatives. This method should fill the _interface_traction and _dinterface_traction_djump varaibles
  virtual void computeInterfaceTractionAndDerivatives() = 0;

  /// Base name of the material system
  const std::string _base_name;

  /// the value of the traction in local coordinates
  MaterialProperty<RealVectorValue> & _interface_traction;

  /// the traction's derivatives wrt the displacement jump in local coordinates
  MaterialProperty<RankTwoTensor> & _dinterface_traction_djump;

  /// the traction's derivatives wrt the fluid velocity jump in local coordinates
  MaterialProperty<RankTwoTensor> & _dinterface_traction_djump_vf;

  /// the traction's derivatives wrt the fluid velocity jump in local coordinates
  MaterialProperty<RealVectorValue> & _dinterface_traction_dpressure;

  /// The displacment jump in local coordaintes
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump;

  /// the value of the traction in local coordinates
  MaterialProperty<Real> & _interface_pressure;

  /// the pressure's derivatives wrt the displacement jump in local coordinates
  MaterialProperty<RealVectorValue> & _dinterface_pressure_djump;

  /// the pressure's derivatives wrt the fluid velocity jump in local coordinates
  MaterialProperty<RealVectorValue> & _dinterface_pressure_djump_vf;

  /// The displacment jump in local coordaintes
  const MaterialProperty<RealVectorValue> & _interface_fluid_vel_jump;

};
