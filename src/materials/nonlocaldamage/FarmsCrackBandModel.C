//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*
Elk Compute Smeared Cracking Stress Model
Created by Chunhui Zhao, Oct 15th, 2024
Rewrite the smeared crack model, add energy regularization
Regularization takes place on equvalent strain

- Pure Solid Mechanics
- Take regularizated equvalent strain as input

*/

#include "FarmsCrackBandModel.h"
#include "ElasticityTensorTools.h"
#include "StressUpdateBase.h"
#include "Conversion.h"

registerMooseObject("farmsApp", FarmsCrackBandModel);

InputParameters
FarmsCrackBandModel::validParams()
{
  InputParameters params = ComputeMultipleInelasticStress::validParams();
  params.addClassDescription("Compute stress using a crack band model");
  params.addRequiredCoupledVar(
      "cracking_stress",
      "The stress threshold beyond which cracking occurs. Negative values prevent cracking.");
  params.addRequiredParam<Real>("hb", "Crack-band size (element length)");
  params.addRequiredParam<Real>("Gf", "Fracture energy G_f");
  params.set<std::vector<MaterialName>>("inelastic_models") = {};
  return params;
}

FarmsCrackBandModel::FarmsCrackBandModel(const InputParameters & parameters)
  : ComputeMultipleInelasticStress(parameters),
  _cracking_stress(coupledValue("cracking_stress")),
  _crack_damage(declareProperty<Real>(_base_name + "crack_damage")),
  _crack_damage_old(getMaterialPropertyOld<Real>(_base_name + "crack_damage")),
  _eqstrain(declareProperty<Real>("eqstrain")),
  _eqstrain_old(getMaterialPropertyOld<Real>("eqstrain")),
  _kappa(declareProperty<Real>("eqstrain_max")),
  _kappa_old(getMaterialPropertyOld<Real>("eqstrain_max")),
  _crack_rotation(declareProperty<RankTwoTensor>(_base_name + "crack_rotation")),
  _crack_rotation_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "crack_rotation")),
  _hb(getParam<Real>("hb")),
  _Gf(getParam<Real>("Gf"))
{
}

void
FarmsCrackBandModel::initQpStatefulProperties()
{
  _crack_damage[_qp] = 0.0;
  _eqstrain[_qp] = 0.0;
  _kappa[_qp] = 0.0;
  _crack_rotation[_qp] = RankTwoTensor::Identity();
}

void
FarmsCrackBandModel::computeQpStress()
{
  bool force_elasticity_rotation = false;

  //(0) get strain at onset of strength criterion
  const Real youngs_modulus =
  ElasticityTensorTools::getIsotropicYoungsModulus(_elasticity_tensor[_qp]);
  Real cracking_strain = _cracking_stress[_qp] / youngs_modulus;

  //(1) update strain tensor
  _elastic_strain[_qp] = _elastic_strain_old[_qp] + _strain_increment[_qp];

  //(2) compute the equivalent strain
  //type 1: Mazars
  RealVectorValue strain_in_crack_dir;
  computeCrackStrainAndOrientation(strain_in_crack_dir);
  Real strain_dir0_positive = std::max(strain_in_crack_dir(0), 0.0);
  Real strain_dir1_positive = std::max(strain_in_crack_dir(1), 0.0);
  Real strain_dir2_positive = std::max(strain_in_crack_dir(2), 0.0);
  Real eqstrain = std::sqrt(strain_dir0_positive*strain_dir0_positive+strain_dir1_positive*strain_dir1_positive+strain_dir2_positive*strain_dir2_positive);
  _eqstrain[_qp] = eqstrain;

  //type 2: Rankine
  // RealVectorValue stress_in_crack_dir;
  // computeCrackStressAndOrientation(stress_in_crack_dir);
  // Real stress_dir0_positive = std::max(stress_in_crack_dir(0), 0.0);
  // Real stress_dir1_positive = std::max(stress_in_crack_dir(1), 0.0);
  // Real stress_dir2_positive = std::max(stress_in_crack_dir(2), 0.0);
  // Real stress_max = std::sqrt(stress_dir0_positive*stress_dir0_positive+stress_dir1_positive*stress_dir1_positive+stress_dir2_positive*stress_dir2_positive);
  // Real eqstrain = stress_max / youngs_modulus;
  // _eqstrain[_qp] = eqstrain;

  //(3) update kappa
  Real k = std::max(_kappa_old[_qp], eqstrain);
  _kappa[_qp] = k;

  //(4) compute ft and mesh-scaled epsf
  Real ft = _cracking_stress[_qp];
  Real epsf = 0.5 * cracking_strain + _Gf/( _hb * ft ); 
  // ensure epsf ≥ eps0 to avoid snap‐back
  epsf = std::max(epsf, cracking_strain);

  //(5) exponential damage law ω(κ) from Eq. (10)
  Real omega = 0.0;
  if (k > cracking_strain)
  {
    Real arg = -(k - cracking_strain)/(epsf - cracking_strain);
    omega = 1.0 - (cracking_strain/k) * std::exp(arg);
  }
  _crack_damage[_qp] = omega;

  //(6) compute damage stress and jacobian
  RankTwoTensor damage_stress = (1 - _crack_damage[_qp]) * _elasticity_tensor[_qp] * _elastic_strain[_qp];
  _stress[_qp] = damage_stress;
  _Jacobian_mult[_qp] = (1 - _crack_damage[_qp]) * _elasticity_tensor[_qp];
  force_elasticity_rotation = true;

  if (_perform_finite_strain_rotations)
  {
    finiteStrainRotation(force_elasticity_rotation);
    _crack_rotation[_qp] = _rotation_increment[_qp] * _crack_rotation[_qp];
  }
}

//compute principal strains and rotation tensor
void
FarmsCrackBandModel::computeCrackStrainAndOrientation(
    RealVectorValue & strain_in_crack_dir)
{
  // The rotation tensor is ordered such that directions for pre-existing cracks appear first
  // in the list of columns.  For example, if there is one existing crack, its direction is in the
  // first column in the rotation tensor.

  std::vector<Real> eigval(3, 0.0);
  RankTwoTensor eigvec;

  _elastic_strain[_qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);

  // If the elastic strain is beyond the cracking strain, save the eigen vectors as
  // the rotation tensor. Reverse their order so that the third principal strain
  // (most tensile) will correspond to the first crack.
  _crack_rotation[_qp].fillColumn(0, eigvec.column(2));
  _crack_rotation[_qp].fillColumn(1, eigvec.column(1));
  _crack_rotation[_qp].fillColumn(2, eigvec.column(0));

  strain_in_crack_dir(0) = eigval[2];
  strain_in_crack_dir(1) = eigval[1];
  strain_in_crack_dir(2) = eigval[0];
}

//compute principal stress in crack direction and rotation tensor
void
FarmsCrackBandModel::computeCrackStressAndOrientation(
    RealVectorValue & stress_in_crack_dir)
{
  // The rotation tensor is ordered such that directions for pre-existing cracks appear first
  // in the list of columns.  For example, if there is one existing crack, its direction is in the
  // first column in the rotation tensor.

  std::vector<Real> eigval(3, 0.0);
  RankTwoTensor eigvec;

  _stress_old[_qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);

  // If the elastic strain is beyond the cracking strain, save the eigen vectors as
  // the rotation tensor. Reverse their order so that the third principal strain
  // (most tensile) will correspond to the first crack.
  _crack_rotation[_qp].fillColumn(0, eigvec.column(2));
  _crack_rotation[_qp].fillColumn(1, eigvec.column(1));
  _crack_rotation[_qp].fillColumn(2, eigvec.column(0));

  stress_in_crack_dir(0) = eigval[2];
  stress_in_crack_dir(1) = eigval[1];
  stress_in_crack_dir(2) = eigval[0];
}