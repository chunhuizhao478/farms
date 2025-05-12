//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "NDSmallDeformationElasticityModel.h"
#include "NDSmallDeformationPlasticityModel.h"

InputParameters
NDSmallDeformationElasticityModel::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();

  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");

  return params;
}

NDSmallDeformationElasticityModel::NDSmallDeformationElasticityModel(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _plasticity_model(nullptr),
    _elastic_strain(declareProperty<RankTwoTensor>(prependBaseName("elastic_strain")))
{
}

void
NDSmallDeformationElasticityModel::setQp(unsigned int qp)
{
  _qp = qp;
  if (_plasticity_model)
    _plasticity_model->setQp(qp);
}

void
NDSmallDeformationElasticityModel::setPlasticityModel(
    NDSmallDeformationPlasticityModel * plasticity_model)
{
  _plasticity_model = plasticity_model;
  _plasticity_model->setElasticityModel(this);
}

void
NDSmallDeformationElasticityModel::initQpStatefulProperties()
{
  _elastic_strain[_qp].zero();
}

void
NDSmallDeformationElasticityModel::updateState(const RankTwoTensor & mechanical_strain,
                                             RankTwoTensor & stress)
{
  _elastic_strain[_qp] = mechanical_strain;

  if (_plasticity_model)
    _plasticity_model->updateState(stress, _elastic_strain[_qp]);
  else
    stress = computeStress(_elastic_strain[_qp]);
}

void
NDSmallDeformationElasticityModel::updateStateDF(const RankTwoTensor & mechanical_strain,
                                                 RankTwoTensor & stress,
                                                 RankFourTensor & Jacobian_mult)
{
  _elastic_strain[_qp] = mechanical_strain;
  
  stress = computeStress(_elastic_strain[_qp]);

  Jacobian_mult = computeJacobian(_elastic_strain[_qp]);

}

