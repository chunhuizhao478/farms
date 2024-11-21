#pragma once

#include "InterfaceKernel.h"
#include "JvarMapInterface.h"

/// Base class for implementing DG cohesive zone models (CZM) for 1D,2D, and 3D
/// traction separation laws. This kernel operates only on
/// a single displacement and fluid velocity compenent.
/// One kernel is required for each displacement component.
class PoroCZMInterfaceKernelPressurebase : public JvarMapKernelInterface<InterfaceKernel>
{
public:
  static InputParameters validParams();
  PoroCZMInterfaceKernelPressurebase(const InputParameters & parameters);

protected:
  Real computeQpResidual(Moose::DGResidualType type) override;
  Real computeQpJacobian(Moose::DGJacobianType type) override;
  Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  /// method computing the derivative of residual[_component] w.r.t pressure
  virtual Real computeDResidualDPressure(const Moose::DGJacobianType & type) const = 0; 
                                                                       

  /// Base name of the material system that this kernel applies to
  const std::string _base_name;

  /// the displacement component this kernel is operating on (0=x, 1=y, 2 =z)
  const unsigned int _component;

  /// Coupled Pressure variable ID
  ///@{
  unsigned int _p_var_num;
  unsigned int _p_var_num_porepressure_neighbor_var;
  ///@}

  // pointer to pressure variable
  MooseVariable & _p_var;

  const VariablePhiValue & _phi;
  const VariablePhiValue & _phi_neighbor;

  // values of the traction and traction derivatives used
  ///@{
  const MaterialProperty<RealVectorValue> & _dtraction_dpressure_global;
  ///@}

  // values of the traction and traction derivatives used
  ///@{
  ///@}

};