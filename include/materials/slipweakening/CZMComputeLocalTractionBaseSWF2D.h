#pragma once

#include "InterfaceMaterial.h"
/**
This is the base function of CZM, In Slip-Weakening Implementation, we declare additional material properties besides interface (traction/jump):
_jump_track_opening : track the jump value due to fault opening and rejoin (after the rejoin, save the current jump, the jump in friction law now becomes: cumulative jump - _jump_track_opening)
_flag_track_opening : flag used to track the opening, 0 is close, 1 is open.
_jump_track_reversal : track the jump value due to reversal of shear stress (after the shear stress change the sign, save the current jump, the jump in friction law now becomes: cumulative jump - _jump_track_reversal)
_T1 : track the current shear stress (for slip reversal case)
_jump_effective : effective jump value (compared with Dc to determine the strength)
**/
class CZMComputeLocalTractionBaseSWF2D : public InterfaceMaterial
{
public:
  static InputParameters validParams();
  CZMComputeLocalTractionBaseSWF2D(const InputParameters & parameters);

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
  MaterialProperty<Real> & _flag_track_opening;
  MaterialProperty<Real> & _flag_track_activecase;

  MaterialProperty<Real> & _jump_track_opening;
  MaterialProperty<Real> & _jump_track_reversal;

  MaterialProperty<Real> & _T1;
  MaterialProperty<Real> & _T2;

  MaterialProperty<Real> & _jump_effective;
  

};