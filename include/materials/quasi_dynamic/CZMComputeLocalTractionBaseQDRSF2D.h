#pragma once

#include "InterfaceMaterial.h"

class CZMComputeLocalTractionBaseQDRSF2D : public InterfaceMaterial
{
public:
  static InputParameters validParams();
  CZMComputeLocalTractionBaseQDRSF2D(const InputParameters & parameters);

protected:
  void initQpStatefulProperties() override;
  void computeQpProperties() override;

  /// Compute the local traction and derivatives. This method should fill the _interface_traction and _dinterface_traction_djump varaibles
  virtual void computeInterfaceTractionAndDerivatives() = 0;

  /// Base name of the material system
  const std::string _base_name;
  Real _V_o;         // Reference slip rate
  Real _f_o;         // Initial friction coefficient
  Real _a;           // Direct effect parameter 
  Real _b;           // State variable evolution parameter
  Real _L; 
  Real _Vini; //initial velocity 
  Real _statevarini; //initial state variable

  //Additional Material Properties
  MaterialProperty<Real> & _sliprate;
  MaterialProperty<Real> & _slip;
  MaterialProperty<Real> & _statevar;
  MaterialProperty<Real> & _statevar_dot;




};