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

#include "FarmsComputeSmearedCrackingStress.h"
#include "ElasticityTensorTools.h"
#include "StressUpdateBase.h"
#include "Conversion.h"

registerMooseObject("farmsApp", FarmsComputeSmearedCrackingStress);

InputParameters
FarmsComputeSmearedCrackingStress::validParams()
{
  InputParameters params = ComputeMultipleInelasticStress::validParams();
  params.addClassDescription("Compute stress using a fixed smeared cracking model");
  params.addRequiredCoupledVar(
      "cracking_stress",
      "The stress threshold beyond which cracking occurs. Negative values prevent cracking.");
  params.addRequiredParam<Real>("paramA", "parameter used in the damage evolution law");
  params.addRequiredParam<Real>("paramB", "parameter used in the damage evolution law");
  params.set<std::vector<MaterialName>>("inelastic_models") = {};
  return params;
}

FarmsComputeSmearedCrackingStress::FarmsComputeSmearedCrackingStress(const InputParameters & parameters)
  : ComputeMultipleInelasticStress(parameters),
    _cracking_stress(coupledValue("cracking_stress")),
    _crack_damage(declareProperty<Real>(_base_name + "crack_damage")),
    _crack_damage_old(getMaterialPropertyOld<Real>(_base_name + "crack_damage")),
    _eqstrain_local(declareProperty<Real>("eqstrain_local")),
    _eqstrain_local_old(getMaterialPropertyOld<Real>("eqstrain_local")),
    _eqstrain_nonlocal(getMaterialProperty<Real>("eqstrain_nonlocal")),
    _kappa(declareProperty<Real>("eqstrain_max")),
    _kappa_old(getMaterialPropertyOld<Real>("eqstrain_max")),
    _crack_rotation(declareProperty<RankTwoTensor>(_base_name + "crack_rotation")),
    _crack_rotation_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "crack_rotation")),
    _paramA(getParam<Real>("paramA")),
    _paramB(getParam<Real>("paramB"))
{
}

void
FarmsComputeSmearedCrackingStress::initQpStatefulProperties()
{
  _crack_damage[_qp] = 0.0;
  _eqstrain_local[_qp] = 0.0;
  _kappa[_qp] = 0.0;
  _crack_rotation[_qp] = RankTwoTensor::Identity();
}

void
FarmsComputeSmearedCrackingStress::computeQpStress()
{
  // (0) Elastic strain update
  _elastic_strain[_qp] = _elastic_strain_old[_qp] + _strain_increment[_qp];

  // (1) Retrieve material parameters and compute cracking strain ε₀
  const Real E    = ElasticityTensorTools::getIsotropicYoungsModulus(_elasticity_tensor[_qp]);
  const Real eps0 = _cracking_stress[_qp] / E;

  // (2) Compute Mazars‐type equivalent strain ε̃ and principal directions
  RealVectorValue eps_dir;
  computeCrackStrainAndOrientation(eps_dir);
  Real eps_dir0 = std::max(eps_dir(0), 0.0);
  Real eps_dir1 = std::max(eps_dir(1), 0.0);
  Real eps_dir2 = std::max(eps_dir(2), 0.0);
  Real eqstrain_local = std::sqrt(eps_dir0*eps_dir0 + eps_dir1*eps_dir1 + eps_dir2*eps_dir2);
  _eqstrain_local[_qp] = eqstrain_local;

  // (3) Update history κ = max(κ_old, ε̃)
  Real kappa = std::max(_kappa_old[_qp], _eqstrain_nonlocal[_qp]);
  _kappa[_qp] = kappa;

  //(4) Exponential damage law ω(κ)
  Real omega = 0.0;
  if (kappa > eps0)
  {
    Real arg1 = std::exp(-_paramB * (kappa - eps0));
    Real arg2 = 1 - _paramA + _paramA * arg1;
    omega = 1.0 - eps0 / kappa * arg2;
  }
  _crack_damage[_qp] = omega;

  // (5) Build consistent tangent and stress
  const RankFourTensor & De  = _elasticity_tensor[_qp];
  const RankTwoTensor  & eps = _elastic_strain[_qp];
  RankTwoTensor De_eps = De * eps;

  // // (6a) dω/dκ from exponential law
  // Real domega_dk = 0.0;
  // if (kappa > eps0 + 1e-12)    // tiny epsilon to avoid exact-zero division
  // {}

  // // (6b) Loading flag: only active when ε̃ > κ_old
  // Real loading = (_eqstrain_nonlocal[_qp] > _kappa_old[_qp] ? 1.0 : 0.0);

  // // (6c) ∂ε̃/∂ε via Mazars
  // RankTwoTensor depsde; 
  // depsde.zero();
  // if (_eqstrain_nonlocal[_qp] > 0.0)
  // {
  //   for (unsigned i = 0; i < 3; ++i)
  //   {
  //     Real eps_i_pos = std::max(eps_dir(i), 0.0);
  //     if (eps_i_pos > 0.0)
  //     {
  //       const RealVectorValue ni = _crack_rotation[_qp].column(i);
  //       depsde += (eps_i_pos/_eqstrain_nonlocal[_qp]) * RankTwoTensor::outerProduct(ni, ni);
  //     }
  //   }
  // }
  // RankTwoTensor domega_de = domega_dk * loading * depsde;

  // usingTensorIndices(i,j,k,l);  // bring i,j,k,l into scope
  // RankFourTensor soft = De_eps.template times<0,1,2,3>(domega_de);

  // (6d) Consistent tangent: (1-ω)De  −  (De:ε) ⊗ (∂ω/∂ε)
  RankFourTensor tangent = (1.0 - omega)*De;

  // (6e) Assign stress and Jacobian multiplier
  _stress[_qp]        = (1.0 - omega)*De_eps;
  _Jacobian_mult[_qp] = tangent;

  // (7) Finite‐strain rotation if needed
  if (_perform_finite_strain_rotations)
  {
    finiteStrainRotation(true);
    _crack_rotation[_qp] = _rotation_increment[_qp] * _crack_rotation[_qp];
  }
}

void
FarmsComputeSmearedCrackingStress::computeCrackStrainAndOrientation(
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
