#pragma once

#include "InterfaceMaterial.h"
/**
This is the base function of CZM, In Slip Weakening Friction, we track the accumulated slip

accumulated_slip_along_normal: accumulated slip along normal direction
accumulated_slip_along_strike: accumulated slip along strike direction
accumulated_slip_along_dip: accumulated slip along dip direction

slip_along_normal: slip along normal
slip_along_strike: slip along strike
slip_along_dip: slip along dip

**/
class CZMComputeLocalTractionBaseLSW3D : public InterfaceMaterial
{
public:
  static InputParameters validParams();
  CZMComputeLocalTractionBaseLSW3D(const InputParameters & parameters);

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

  /// The displacment jump in local coordaintes
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump;

  //Additional Material Properties
  MaterialProperty<Real> & _accumulated_slip_along_normal;
  MaterialProperty<Real> & _accumulated_slip_along_strike;
  MaterialProperty<Real> & _accumulated_slip_along_dip;

  MaterialProperty<Real> & _slip_along_normal;
  MaterialProperty<Real> & _slip_along_strike;
  MaterialProperty<Real> & _slip_along_dip;

  MaterialProperty<Real> & _jump_track_dip;
  MaterialProperty<Real> & _T3;

};