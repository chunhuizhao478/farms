//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Material.h"
#include "RankTwoTensorForward.h"
#include "BaseNameInterface.h"

class NDSmallDeformationPlasticityModel;

class NDSmallDeformationElasticityModel : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  NDSmallDeformationElasticityModel(const InputParameters & parameters);

  /// Set the current quadrature point
  virtual void setQp(unsigned int qp);

  /// Set the associated plasticity model
  virtual void setPlasticityModel(NDSmallDeformationPlasticityModel * plasticity_model);

  /**
   * Compute the stress given the mechanical strain. Also performs the plasticity update, if any.
   * @param mechanical_strain The mechanical strain
   * @param stress            The stress
   */
  virtual void updateState(const RankTwoTensor & mechanical_strain, RankTwoTensor & stress);

  //use non-AD implementation
  virtual void updateStateDF(const RankTwoTensor & mechanical_strain,
                             RankTwoTensor & stress,
                             RankFourTensor & Jacobian_mult);

  /**
   * The stress-strain relation
   * @param strain The given strain
   * @return The computed stress given the strain and the constitutive relation
   */
  virtual RankTwoTensor computeStress(const RankTwoTensor & strain) = 0;
  
  //use non-AD implementation
  //compute jacobian of the stress w/r/t strain
  virtual RankFourTensor computeJacobian(const RankTwoTensor & strain) = 0;

  // @{ Retained as empty methods to avoid a warning from Material.C in framework. These methods are
  // unused in all inheriting classes and should not be overwritten.
  void resetQpProperties() final {}
  void resetProperties() final {}
  // @}

  //helper function, grab from RaccoonUtils, make it non-AD
  Real Macaulay(const Real x, const bool deriv = false);

  std::vector<Real> Macaulay(const std::vector<Real> & v, const bool deriv = false);

  RankTwoTensor spectralDecomposition(const RankTwoTensor & r2t);

protected:
  virtual void initQpStatefulProperties() override;

  /// The optional plasticity model
  NDSmallDeformationPlasticityModel * _plasticity_model;

  /// The elastic strain
  MaterialProperty<RankTwoTensor> & _elastic_strain;
};
