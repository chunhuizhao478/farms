//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Material.h"
#include "RankTwoTensorForward.h"
#include "SingleVariableReturnMappingSolution.h"
#include "BaseNameInterface.h"
#include "PlasticHardeningModel.h"

class NDSmallDeformationElasticityModel;

class NDSmallDeformationPlasticityModel : public Material,
                                        public SingleVariableReturnMappingSolution,
                                        public BaseNameInterface
{
public:
  static InputParameters validParams();

  NDSmallDeformationPlasticityModel(const InputParameters & parameters);

  virtual void initialSetup() override;

  /// Set the current quadrature point
  virtual void setQp(unsigned int qp);

  /// Set the associated elasticity model
  virtual void setElasticityModel(NDSmallDeformationElasticityModel * elasticity_model);

  /**
   * Update the stress and elastic strain if need to following the specified plastic flow
   * @param stress         The stress
   * @param elastic_strain The elastic strain
   */
  virtual void updateState(RankTwoTensor & stress, RankTwoTensor & elastic_strain) = 0;

  // @{ Retained as empty methods to avoid a warning from Material.C in framework. These methods are
  // unused in all inheriting classes and should not be overwritten.
  void resetQpProperties() final {}
  void resetProperties() final {}
  // @}

protected:
  virtual void initQpStatefulProperties() override;

  /// The elasticity model
  NDSmallDeformationElasticityModel * _elasticity_model;

  /// The plastic strain
  MaterialProperty<RankTwoTensor> & _plastic_strain;
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;

  /// The (scalar) effective plastic strain
  MaterialProperty<Real> & _ep;
  const MaterialProperty<Real> & _ep_old;

  /// The flow direction
  MaterialProperty<RankTwoTensor> & _Np;

  /// The plastic hardening model
  PlasticHardeningModel * _hardening_model;
};
