
#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"

/**
 * DarcyFluidVelocity 
 */
class DarcyFluidVelocity : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  DarcyFluidVelocity(const InputParameters & parameters);

  void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void VelocityIntegrityCheck();

  /// Coupled Velocities variables
  unsigned int _nvel;

  /// Velocities variables
  std::vector<const VariableValue *> _vel;

  /// Gradient of Velocities
  std::vector<const VariableGradient *> _grad_vel;

    /// Base name of the material system
  const std::string _base_name;

  MaterialProperty<RankTwoTensor> & _darcy_vel_grad;

  const Real & _current_elem_volume;
};