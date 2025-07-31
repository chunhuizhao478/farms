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

#include "FarmsComputeSmearedCrackingStressGrads.h"
#include "ElasticityTensorTools.h"
#include "StressUpdateBase.h"
#include "Conversion.h"

registerMooseObject("farmsApp", FarmsComputeSmearedCrackingStressGrads);

InputParameters
FarmsComputeSmearedCrackingStressGrads::validParams()
{
  InputParameters params = ComputeMultipleInelasticStress::validParams();
  params.addClassDescription("Compute stress using a fixed smeared cracking model");
  params.addRequiredCoupledVar(
      "cracking_stress",
      "The stress threshold beyond which cracking occurs. Negative values prevent cracking.");
  params.addRequiredCoupledVar(
      "nonlocal_eqstrain",
      "The nonlocal equivalent strain used in the damage evolution law");
  params.addRequiredParam<Real>("paramA", "parameter used in the damage evolution law");
  params.addRequiredParam<Real>("paramB", "parameter used in the damage evolution law");
  params.set<std::vector<MaterialName>>("inelastic_models") = {};

  //add initial damage
  params.addRequiredCoupledVar(
      "initial_crack_damage",
      "Initial damage for crack_damage material property");
  
  params.addParam<bool>("porous_flow_coupling", false, "Enable porous flow coupling");
  params.addParam<Real>("intrinsic_permeability", 5e-19, "Intrinsic permeability in m^2");
  //Permeability models
  //Exponential permeability model
  params.addParam<bool>("exponential_permeability_model", false,
                        "Use an exponential function for the effective permeability");
  params.addParam<Real>("coeff_b", -1.0,
                        "Coefficient for the exponential function in the effective permeability");
  //Darcy-Poiseuille permeability model
  params.addParam<bool>("darcy_poiseuille_permeability_model",
                        false,
                        "Use Darcy-Poiseuille model for the effective permeability");
  params.addParam<Real>("wc", -1.0, "ultimate crack width for Darcy-Poiseuille model");
  params.addParam<Real>("perm_exponent", -1.0,
                        "Exponent for the Darcy-Poiseuille model for the effective permeability");
  return params;
}

FarmsComputeSmearedCrackingStressGrads::FarmsComputeSmearedCrackingStressGrads(const InputParameters & parameters)
  : ComputeMultipleInelasticStress(parameters),
    _cracking_stress(coupledValue("cracking_stress")),
    _crack_damage(declareProperty<Real>(_base_name + "crack_damage")),
    _crack_damage_old(getMaterialPropertyOld<Real>(_base_name + "crack_damage")),
    _eqstrain_local(declareProperty<Real>("eqstrain_local")),
    _eqstrain_local_old(getMaterialPropertyOld<Real>("eqstrain_local")),
    _eqstrain_nonlocal(coupledValue("nonlocal_eqstrain")),
    _kappa(declareProperty<Real>("eqstrain_max")),
    _kappa_old(getMaterialPropertyOld<Real>("eqstrain_max")),
    _crack_rotation(declareProperty<RankTwoTensor>(_base_name + "crack_rotation")),
    _crack_rotation_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "crack_rotation")),
    _paramA(getParam<Real>("paramA")),
    _paramB(getParam<Real>("paramB")),
    _initial_crack_damage(coupledValue("initial_crack_damage")),
    //porous flow coupling
    _porous_flow_coupling(getParam<bool>("porous_flow_coupling")),
    _intrinsic_permeability(getParam<Real>("intrinsic_permeability")),
    // define effective permeability
    _effective_perm(declareProperty<RealTensorValue>("effective_perm")),
    _effective_perm_old(getMaterialPropertyOldByName<RealTensorValue>("effective_perm")),
    // Exponential permeability model
    _exponential_permeability_model(getParam<bool>("exponential_permeability_model")),
    _coeff_b(getParam<Real>("coeff_b")),
    // Darcy-Poiseuille permeability model
    _darcy_poiseuille_permeability_model(getParam<bool>("darcy_poiseuille_permeability_model")),
    _wc(getParam<Real>("wc")),
    _perm_exponent(getParam<Real>("perm_exponent"))
{
}

void
FarmsComputeSmearedCrackingStressGrads::initQpStatefulProperties()
{
  _crack_damage[_qp] = _initial_crack_damage[_qp];
  _eqstrain_local[_qp] = 0.0;
  _kappa[_qp] = 0.0;
  _crack_rotation[_qp] = RankTwoTensor::Identity();
}

void
FarmsComputeSmearedCrackingStressGrads::computeQpStress()
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
    // Real arg1 = std::exp(-_paramB * (kappa - eps0));
    // Real arg2 = 1 - _paramA + _paramA * arg1;
    // omega = 1.0 - eps0 / kappa * arg2;

    omega = 1.0 - eps0 / kappa * (1.0 - _paramA) - _paramA / std::exp(_paramB * (kappa - eps0));

  }

  // irreversible crack damage, set to initial damage if it is smaller
  if (omega < _initial_crack_damage[_qp])
  {
    omega = _initial_crack_damage[_qp];
  }

  // Ensure damage is non-decreasing (enforce irreversibility)
  if (omega < _crack_damage_old[_qp])
  {
    omega = _crack_damage_old[_qp];
  }

  //save the damage
  _crack_damage[_qp] = omega;

  // (5) Build consistent tangent and stress
  const RankFourTensor & De  = _elasticity_tensor[_qp];
  const RankTwoTensor  & eps = _elastic_strain[_qp];
  RankTwoTensor De_eps = De * eps;

  // (6d) Consistent tangent: (1-ω)De - (De:ε) ⊗ (∂ω/∂ε)
  RankFourTensor tangent = (1.0 - omega + 0.01) * De;

  // (6e) Assign stress and Jacobian multiplier
  _stress[_qp] = (1.0 - omega + 0.01) * De_eps;
  _Jacobian_mult[_qp] = tangent;

  // (7) Finite‐strain rotation if needed
  if (_perform_finite_strain_rotations)
  {
    finiteStrainRotation(true);
    _crack_rotation[_qp] = _rotation_increment[_qp] * _crack_rotation[_qp];
  }

  // Compute effective permeability
  updatePermeabilityForCracking();
}

void
FarmsComputeSmearedCrackingStressGrads::computeCrackStrainAndOrientation(
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

void
FarmsComputeSmearedCrackingStressGrads::updatePermeabilityForCracking()
{

  // If porous flow coupling is not enabled, return
  if (!_porous_flow_coupling)
    return;

  // Get transformation matrix
  const RankTwoTensor & R = _crack_rotation[_qp];

  // Initialize effective permeability new
  RankTwoTensor effective_perm_new;

  //Compute the intrinsic permeability
  RankTwoTensor perm_intrinsic = _intrinsic_permeability * RankTwoTensor::Identity();

  // Initialize effective permeability new
  // exponential permeability model 
  if (_exponential_permeability_model){
    effective_perm_new = perm_intrinsic * std::exp( _crack_damage[_qp] * _coeff_b );
  }
  // darcy-poiseuille permeability model
  else if (_darcy_poiseuille_permeability_model){
    //Compute crack opening
    //wc is the ultimate crack opening
    Real w = _crack_damage[_qp] * _wc; 

    //Compute permeability in the damage zone
    RankTwoTensor kf = std::pow(w, 2) / (12.0) * RankTwoTensor::Identity();

    //Compute permeability
    effective_perm_new = perm_intrinsic + std::pow(_crack_damage[_qp], _perm_exponent) * (kf - perm_intrinsic);
  }
  else {
    mooseError("Unknown permeability model type.");
  }

  // Rotate back to global frame
  effective_perm_new.rotate(R);

  // Update effective perm
  _effective_perm[_qp] = effective_perm_new;

}
