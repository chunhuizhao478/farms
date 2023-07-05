#pragma once

#include "InterfaceMaterial.h"
/**
This is the base function of CZM, In Rate-and-State Implementation, we declare additional material properties besides interface (traction/jump):
_alongfaultvel_strike_plus   : velocity along the upper fault surface in strike direction
_alongfaultvel_strike_minus  : velocity along the lower fault surface in strike direction
_alongfaultvel_dip_plus      : velocity along the upper fault surface in dip direction
_alongfaultvel_dip_minus     : velocity along the lower fault surface in dip direction
_alongfaultvel_normal_plus   : velocity along the upper fault surface in normal direction
_alongfaultvel_normal_minus  : velocity along the lower fault surface in normal direction

_alongfaultdisp_strike_plus  : displacement along the upper fault surface in strike direction
_alongfaultdisp_strike_minus : displacement along the lower fault surface in strike direction
_alongfaultdisp_dip_plus     : displacement along the upper fault surface in dip direction
_alongfaultdisp_dip_minus    : displacement along the lower fault surface in dip direction
_alongfaultdisp_normal_plus  : displacement along the upper fault surface in normal direction
_alongfaultdisp_normal_minus : displacement along the lower fault surface in normal direction

_sliprate_strike : slip rate along the fault in strike direction
_sliprate_dip    : slip rate along the fault in dip direction
_sliprate_mag    : slip rate magnitude along the fault combining both strike and dip direction

_statevar : state variable 

_traction_strike : traction in the strike direction
_traction_dip    : traction in the dip direction
_traction_normal : traction in the normal direction
**/
class CZMComputeLocalTractionBaseRSF2D : public InterfaceMaterial
{
public:
  static InputParameters validParams();
  CZMComputeLocalTractionBaseRSF2D(const InputParameters & parameters);

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
  MaterialProperty<Real> & _alongfaultvel_strike_plus;
  MaterialProperty<Real> & _alongfaultvel_strike_minus;

  MaterialProperty<Real> & _alongfaultvel_normal_plus;
  MaterialProperty<Real> & _alongfaultvel_normal_minus;
  
  MaterialProperty<Real> & _alongfaultdisp_strike_plus;
  MaterialProperty<Real> & _alongfaultdisp_strike_minus;

  MaterialProperty<Real> & _alongfaultdisp_normal_plus;
  MaterialProperty<Real> & _alongfaultdisp_normal_minus;

  MaterialProperty<Real> & _sliprate_strike;
  MaterialProperty<Real> & _slip_strike;

  MaterialProperty<Real> & _sliprate_normal;
  MaterialProperty<Real> & _slip_normal;

  MaterialProperty<Real> & _sliprate_mag;
  MaterialProperty<Real> & _slip_mag;

  MaterialProperty<Real> & _statevar;

  MaterialProperty<Real> & _traction_strike;
  MaterialProperty<Real> & _traction_normal;

  MaterialProperty<Real> & _sliprate_predict;
  MaterialProperty<Real> & _slip_predict;

  Real _Ts_o; //initial shear traction
  Real _Tn_o; //initial normal traction
  Real _Vini; //initial velocity 
  Real _statevarini; //initial state variable

  //
  MaterialProperty<Real> & _alongfaultvel_x_plus;
  MaterialProperty<Real> & _alongfaultvel_x_minus;

  MaterialProperty<Real> & _alongfaultvel_y_plus;
  MaterialProperty<Real> & _alongfaultvel_y_minus;
  
  MaterialProperty<Real> & _alongfaultdisp_x_plus;
  MaterialProperty<Real> & _alongfaultdisp_x_minus;

  MaterialProperty<Real> & _alongfaultdisp_y_plus;
  MaterialProperty<Real> & _alongfaultdisp_y_minus;

};