#pragma once

#include "InterfaceMaterial.h"

#define usingCZMComputeDarcyFluidVelocityJumpBaseMembers                                     \
  usingInterfaceMaterialMembers;                                                             \
  using CZMComputeDarcyFluidVelocityJumpBase::_base_name;                                    \
  using CZMComputeDarcyFluidVelocityJumpBase::_nvel;                                         \
  using CZMComputeDarcyFluidVelocityJumpBase::_fluid_vel;                                    \
  using CZMComputeDarcyFluidVelocityJumpBase::_fluid_vel_neighbor;                           \
  using CZMComputeDarcyFluidVelocityJumpBase::_fluid_vel_jump_global;                        \
  using CZMComputeDarcyFluidVelocityJumpBase::_interface_fluid_vel_jump;                     \
  using CZMComputeDarcyFluidVelocityJumpBase::_czm_total_rotation

/**
 * This interface material class computes the fluid velocity jump in the interface natural coordinate
 * system. The transformation between local and global coordinates shall be defined in
 * computeLocalDisplacementJump.
 */

class CZMComputeDarcyFluidVelocityJumpBase : public InterfaceMaterial
{
public:
  static InputParameters validParams();
  CZMComputeDarcyFluidVelocityJumpBase(const InputParameters & parameters);

protected:
  void computeQpProperties() override;
  void initQpStatefulProperties() override;

  /// method used to compute the disaplcement jump in interface coordinates according to a
  ///  specific kinematic formulation
  virtual void computeLocalFluidVelocityJump() = 0;

  /// method computing the required rotation matrices
  virtual void computeRotationMatrices();

  /// Base name of the material system
  const std::string _base_name;

  /// number of displacement components
  const unsigned int _nvel;

  //   /// interface type
  // enum class InterfaceType
  // {
  //   Permeable,
  //   ImPermeable,
  //   DynamicDarcy,
  // } _interface_type;

  /// the coupled displacement and neighbor displacement values
  ///@{
  std::vector<const GenericVariableValue<false> *> _fluid_vel;
  std::vector<const GenericVariableValue<false> *> _fluid_vel_neighbor;
  ///@}

  /// the displacement jump in global and interface coordiantes
  ///@{
  GenericMaterialProperty<RealVectorValue, false> & _fluid_vel_jump_global;
  GenericMaterialProperty<RealVectorValue, false> & _interface_fluid_vel_jump;
  ///@}

  // enum class InterfaceType
  // {
  //   Permeable,
  //   Impermeable,
  //   Darcy,
  //   Hydrofault
  // } _interface_type;

  /// the rotation matrix transforming from the interface to the global coordinate systems
  GenericMaterialProperty<RankTwoTensor, false> & _czm_total_rotation;
};