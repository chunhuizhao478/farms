//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "NDSmallDeformationPlasticityModel.h"
#include "NDSmallDeformationElasticityModel.h"

InputParameters
NDSmallDeformationPlasticityModel::validParams()
{
  InputParameters params = Material::validParams();
  params += SingleVariableReturnMappingSolution::validParams();
  params += BaseNameInterface::validParams();

  params.addRequiredParam<MaterialName>("hardening_model", "Name of the plastic hardening model");

  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");

  return params;
}

NDSmallDeformationPlasticityModel::NDSmallDeformationPlasticityModel(const InputParameters & parameters)
  : Material(parameters),
    SingleVariableReturnMappingSolution(parameters),
    BaseNameInterface(parameters),
    _plastic_strain(declareProperty<RankTwoTensor>(prependBaseName("plastic_strain"))),
    _plastic_strain_old(
        getMaterialPropertyOldByName<RankTwoTensor>(prependBaseName("plastic_strain"))),
    _ep(declareProperty<Real>(prependBaseName("effective_plastic_strain"))),
    _ep_old(getMaterialPropertyOldByName<Real>(prependBaseName("effective_plastic_strain"))),
    _Np(declareProperty<RankTwoTensor>(prependBaseName("flow_direction")))
{
}

void
NDSmallDeformationPlasticityModel::initialSetup()
{
  _hardening_model = dynamic_cast<PlasticHardeningModel *>(&getMaterial("hardening_model"));
  if (!_hardening_model)
    paramError("hardening_model",
               "Plastic hardening model " + getParam<MaterialName>("hardening_model") +
                   " is not compatible with " + name());
}

void
NDSmallDeformationPlasticityModel::setQp(unsigned int qp)
{
  _qp = qp;
  _hardening_model->setQp(qp);
}

void
NDSmallDeformationPlasticityModel::setElasticityModel(
    NDSmallDeformationElasticityModel * elasticity_model)
{
  _elasticity_model = elasticity_model;
}

void
NDSmallDeformationPlasticityModel::initQpStatefulProperties()
{
  _plastic_strain[_qp].zero();
  _ep[_qp] = 0;
}
