#pragma once

#include "InterfaceKernel.h"
#include "JvarMapInterface.h"

/// Base class for implementing DG cohesive zone models (CZM) for 1D,2D, and 3D
/// traction separation laws. This kernel operates only on
/// a single displacement and fluid velocity compenent.
/// One kernel is required for each displacement component.
class PoroCZMInterfaceKernelbase : public JvarMapKernelInterface<InterfaceKernel>
{
public:
  static InputParameters validParams();
  PoroCZMInterfaceKernelbase(const InputParameters & parameters);

protected:
  Real computeQpResidual(Moose::DGResidualType type) override;
  Real computeQpJacobian(Moose::DGJacobianType type) override;
  Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  /// method computing the derivative of residual[_component] w.r.t displacement[component_j]
  virtual Real computeDResidualDDisplacement(const unsigned int & component_j,
                                             const Moose::DGJacobianType & type) const = 0;

  // /// method computing the derivative of residual[_component] w.r.t fluid velocity[component_j]
  // virtual Real computeDResidualDFlux(const unsigned int & component_j,
  //                                    const Moose::DGJacobianType & type) const = 0;

  // /// method computing the derivative of residual[_component] w.r.t pressure
  // virtual Real computeDResidualDPressure(const Moose::DGJacobianType & type) const = 0; 
                                                                       

  /// Base name of the material system that this kernel applies to
  const std::string _base_name;

  /// the displacement component this kernel is operating on (0=x, 1=y, 2 =z)
  const unsigned int _component;

  /// number of displacement components
  const unsigned int _ndisp;

  /// Coupled displacement component variable IDs
  ///@{
  std::vector<unsigned int> _disp_var;
  std::vector<unsigned int> _disp_neighbor_var;
  ///@}

  // pointer to displacement variables
  std::vector<MooseVariable *> _vars;

  // /// number of Fluid velocity components
  // const unsigned int _nvel;

  // /// Coupled Fluid velocity component variable IDs
  // ///@{
  // std::vector<unsigned int> _fluid_vel_var;
  // std::vector<unsigned int> _fluid_vel_neighbor_var;
  // ///@}

  // // pointer to Fluid velocity variables
  // std::vector<MooseVariable *> _fluid_vars;

  // /// Coupled Pressure variable ID
  // ///@{
  // unsigned int _pp_var;
  // unsigned int _pp_neighbor_var;
  // ///@}

  // // pointer to Fluid velocity variables
  // MooseVariable * _pp_vars;


  // values of the traction and traction derivatives used
  ///@{
  const MaterialProperty<RealVectorValue> & _traction_global;
  const MaterialProperty<RankTwoTensor> & _dtraction_djump_global;
  // const MaterialProperty<RankTwoTensor> & _dtraction_djump_global_vf;
  // const MaterialProperty<RealVectorValue> & _dtraction_dpressure_global;
  ///@}

  // values of the traction and traction derivatives used
  ///@{
  const MaterialProperty<Real> & _pressure_global;
  const MaterialProperty<RealVectorValue> & _dpressure_djump_global;
  // const MaterialProperty<RealVectorValue> & _dpressure_djump_global_vf;
  ///@}

};