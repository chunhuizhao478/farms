//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ComputeSmallDeformationStress.h"
#include "SmallDeformationElasticityModel.h"
#include "SmallDeformationPlasticityModel.h"

registerMooseObject("farmsApp", ComputeSmallDeformationStress);

InputParameters
ComputeSmallDeformationStress::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();
  params.addClassDescription("The stress calculator given an elasticity model and a plasticity "
                             "model. Small deformation is assumed. Also computes strain increment.");

  params.addRequiredParam<MaterialName>("elasticity_model",
                                        "Name of the elastic stress-strain constitutive model");
  params.addParam<MaterialName>("plasticity_model", "Name of the plasticity model");
  params.addParam<bool>("compute_strain_increment",
                        true,
                        "Whether to compute strain increment (requires elastic_strain to be "
                        "stateful). Set to false for backward compatibility with old checkpoints "
                        "that don't have elastic_strain as a stateful property.");

  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

ComputeSmallDeformationStress::ComputeSmallDeformationStress(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _mechanical_strain(getADMaterialProperty<RankTwoTensor>(prependBaseName("mechanical_strain"))),
    _elastic_strain(getADMaterialProperty<RankTwoTensor>(prependBaseName("elastic_strain"))),
    _compute_strain_increment(getParam<bool>("compute_strain_increment")),
    _elastic_strain_old(_compute_strain_increment
                            ? &getMaterialPropertyOld<RankTwoTensor>(prependBaseName("elastic_strain"))
                            : nullptr),
    _stress(declareADProperty<RankTwoTensor>(prependBaseName("stress"))),
    _strain_increment(_compute_strain_increment
                          ? &declareADProperty<RankTwoTensor>(prependBaseName("strain_increment"))
                          : nullptr)
{
  if (getParam<bool>("use_displaced_mesh"))
    mooseError("The stress calculator needs to run on the undisplaced mesh.");
}

void
ComputeSmallDeformationStress::initialSetup()
{
  _elasticity_model =
      dynamic_cast<SmallDeformationElasticityModel *>(&getMaterial("elasticity_model"));
  if (!_elasticity_model)
    paramError("elasticity_model",
               "Elasticity model " + getParam<MaterialName>("elasticity_model") +
                   " is not compatible with ComputeSmallDeformationStress");

  _plasticity_model =
      isParamValid("plasticity_model")
          ? dynamic_cast<SmallDeformationPlasticityModel *>(&getMaterial("plasticity_model"))
          : nullptr;
  if (_plasticity_model)
    _elasticity_model->setPlasticityModel(_plasticity_model);
}

void
ComputeSmallDeformationStress::initQpStatefulProperties()
{
  _stress[_qp].zero();
  if (_strain_increment)
    (*_strain_increment)[_qp].zero();
}

void
ComputeSmallDeformationStress::computeQpProperties()
{
  _elasticity_model->setQp(_qp);
  _elasticity_model->updateState(_mechanical_strain[_qp], _stress[_qp]);

  // Compute elastic strain increment for this step: current minus old
  // Note: _elastic_strain is AD, _elastic_strain_old is non-AD (from previous timestep)
  if (_compute_strain_increment && _strain_increment && _elastic_strain_old)
  {
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        (*_strain_increment)[_qp](i, j) =
            _elastic_strain[_qp](i, j) - (*_elastic_strain_old)[_qp](i, j);
  }
}
