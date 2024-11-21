#pragma once

#include "StressDivergenceTensors.h"

/**
 * ViscoelasticStressKernel computes the viscoelastic stress contribution using
 * Kelvin-Voigt model with η_kv/G * stress_rate
 */
class ViscoelasticStressKernel : public StressDivergenceTensors
{
public:
  static InputParameters validParams();
  ViscoelasticStressKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Previous timestep stress
  const MaterialProperty<RankTwoTensor> & _stress_older;
  
  /// Current timestep stress
  const MaterialProperty<RankTwoTensor> & _stress;

  /// Kelvin-Voigt viscosity parameter
  const Real _eta_kv;
  
  /// Shear modulus
  const Real _G;

  /// Dimension of the mesh
  const unsigned int _dim;
};

