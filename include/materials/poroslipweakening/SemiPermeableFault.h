/* 
Material Description of Slip Weakening Friction v3
Generalize the computation of sticking traction using consistent displacement jump and nodal reaction forces
*/

#pragma once

#include "InterfaceMaterial.h"

class SemiPermeableFault : public InterfaceMaterial
{
public:
  static InputParameters validParams();
  SemiPermeableFault(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// Base name of the material system
  const std::string _base_name;

  const MaterialProperty<Real> &_Transmissibility;

  const VariableValue & _interface_pressure_secondary; 
  const VariableValue & _interface_pressure_main;

  MaterialProperty<Real> & _across_fault_flux;

};