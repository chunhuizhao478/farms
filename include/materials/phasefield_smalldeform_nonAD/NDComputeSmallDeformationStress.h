//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Material.h"
#include "RankTwoTensorForward.h"
#include "BaseNameInterface.h"

/*
This class is in the same position as ComputeGeneralStressBase.h
Here we needs to define stress and jacobian as material properties
*/

class NDSmallDeformationElasticityModel;
class NDSmallDeformationPlasticityModel;

/**
 * NDComputeSmallDeformationStress computes the stress under small-strain assumptions
 */
class NDComputeSmallDeformationStress : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  NDComputeSmallDeformationStress(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// The elasticity model
  NDSmallDeformationElasticityModel * _elasticity_model;

  /// The plasticity model
  NDSmallDeformationPlasticityModel * _plasticity_model;

  /// The mechanical strain excluding eigen strains from the total strain
  const MaterialProperty<RankTwoTensor> & _mechanical_strain;

  /// The stress
  MaterialProperty<RankTwoTensor> & _stress;

  /// derivative of stress w.r.t. strain (_dstress_dstrain)
  MaterialProperty<RankFourTensor> & _Jacobian_mult;
};
