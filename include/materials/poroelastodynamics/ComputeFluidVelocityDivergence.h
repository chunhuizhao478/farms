
#pragma once

#include "PorousFlowMaterialVectorBase.h"
#include "RankTwoTensor.h"

/**
 * PorousFlowVolumetricStrain computes volumetric strains, and derivatives thereof
 */
class ComputeFluidVelocityDivergence : public PorousFlowMaterialVectorBase
{
public:
  static InputParameters validParams();

  ComputeFluidVelocityDivergence(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// base name used in the Tensor Mechanics strain calculator
  const std::string _base_name;

  /// Darcy Velocity Gradient 
  const MaterialProperty<RankTwoTensor> & _darcy_vel_grad;

  /// Number of Velocities supplied (1 in 1D, 2 in 2D, 3 in 3D)
  const unsigned int _nvel;

  /// MOOSE variable number of the Darcy Velocity provided
  std::vector<unsigned int> _vel_var_num;

  /// The Darcy Velocity Divergence at the quadpoints
  MaterialProperty<Real> & _darcy_vel_div_qp;

  /**
   * The derivative of the Darcy Velocity Divergence with respect to the porous flow variables.
   * Since the Darcy Velocity Divergence depends on derivatives of porous flow variables.
   * this should be multiplied by _grad_phi in kernels
   */
  MaterialProperty<std::vector<RealGradient>> & _ddarcy_vel_div_qp_dvar;
};