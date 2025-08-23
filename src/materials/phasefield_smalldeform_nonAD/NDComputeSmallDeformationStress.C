//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "NDComputeSmallDeformationStress.h"
#include "NDSmallDeformationElasticityModel.h"
#include "NDSmallDeformationPlasticityModel.h"

registerMooseObject("farmsApp", NDComputeSmallDeformationStress);

InputParameters
NDComputeSmallDeformationStress::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();
  params.addClassDescription("The stress calculator given an elasticity model and a plasticity "
                             "model. Small deformation is assumed.");

  params.addRequiredParam<MaterialName>("elasticity_model",
                                        "Name of the elastic stress-strain constitutive model");
  params.addParam<MaterialName>("plasticity_model", "Name of the plasticity model");

  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

NDComputeSmallDeformationStress::NDComputeSmallDeformationStress(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>(prependBaseName("mechanical_strain"))),
  _elastic_strain(getMaterialProperty<RankTwoTensor>(prependBaseName("elastic_strain"))),
  _elastic_strain_old(getMaterialPropertyOld<RankTwoTensor>(prependBaseName("elastic_strain"))),
    _stress(declareProperty<RankTwoTensor>(prependBaseName("stress"))),
  _Jacobian_mult(declareProperty<RankFourTensor>(prependBaseName("Jacobian_mult"))),
  _strain_increment(declareProperty<RankTwoTensor>(prependBaseName("strain_increment")))
{
  if (getParam<bool>("use_displaced_mesh"))
    mooseError("The stress calculator needs to run on the undisplaced mesh.");
}

void
NDComputeSmallDeformationStress::initialSetup()
{
  _elasticity_model =
      dynamic_cast<NDSmallDeformationElasticityModel *>(&getMaterial("elasticity_model"));
  if (!_elasticity_model)
    paramError("elasticity_model",
               "Elasticity model " + getParam<MaterialName>("elasticity_model") +
                   " is not compatible with NDComputeSmallDeformationStress");

  _plasticity_model =
      isParamValid("plasticity_model")
          ? dynamic_cast<NDSmallDeformationPlasticityModel *>(&getMaterial("plasticity_model"))
          : nullptr;
  if (_plasticity_model)
    _elasticity_model->setPlasticityModel(_plasticity_model);
}

void
NDComputeSmallDeformationStress::initQpStatefulProperties()
{
  _stress[_qp].zero();
  _strain_increment[_qp].zero();
}

void
NDComputeSmallDeformationStress::computeQpProperties()
{
  _elasticity_model->setQp(_qp);
  _elasticity_model->updateStateDF(_mechanical_strain[_qp], 
                                   _stress[_qp],
                                   _Jacobian_mult[_qp]);

  // Compute elastic strain increment for this step: current minus old
  _strain_increment[_qp] = _elastic_strain[_qp] - _elastic_strain_old[_qp];
}
