//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*
AD Farms Compute Smeared Cracking Stress Model
Created by Chunhui Zhao, Oct 15th, 2024
Rewrite the smeared crack model, add energy regularization
Regularization takes place on equvalent strain
Uses automatic differentiation

- Pure Solid Mechanics
- Take regularizated equvalent strain as input
*/

#include "ADFarmsComputeSmearedCrackingStressGrads.h"
#include "ElasticityTensorTools.h"
#include "StressUpdateBase.h"
#include "Conversion.h"

registerADMooseObject("farmsApp", ADFarmsComputeSmearedCrackingStressGrads);

InputParameters
ADFarmsComputeSmearedCrackingStressGrads::validParams()
{
  InputParameters params = ADComputeMultipleInelasticStress::validParams();
  params.addClassDescription("Compute stress using a fixed smeared cracking model with automatic differentiation");
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

ADFarmsComputeSmearedCrackingStressGrads::ADFarmsComputeSmearedCrackingStressGrads(const InputParameters & parameters)
  : ADComputeMultipleInelasticStress(parameters),
    _cracking_stress(adCoupledValue("cracking_stress")),
    _crack_damage(declareADProperty<Real>(_base_name + "crack_damage")),
    _crack_damage_old(getMaterialPropertyOld<Real>(_base_name + "crack_damage")),
    _eqstrain_local(declareADProperty<Real>("eqstrain_local")),
    _eqstrain_local_old(getMaterialPropertyOld<Real>("eqstrain_local")),
    _eqstrain_nonlocal(adCoupledValue("nonlocal_eqstrain")),
    _kappa(declareADProperty<Real>("eqstrain_max")),
    _kappa_old(getMaterialPropertyOld<Real>("eqstrain_max")),
    _crack_rotation(declareADProperty<RankTwoTensor>(_base_name + "crack_rotation")),
    _crack_rotation_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "crack_rotation")),
    _paramA(getParam<Real>("paramA")),
    _paramB(getParam<Real>("paramB")),
    _initial_crack_damage(adCoupledValue("initial_crack_damage")),
    //porous flow coupling
    _porous_flow_coupling(getParam<bool>("porous_flow_coupling")),
    _intrinsic_permeability(getParam<Real>("intrinsic_permeability")),
    // define effective permeability
    _effective_perm(declareADProperty<RealTensorValue>("effective_perm")),
    _effective_perm_old(getMaterialPropertyOldByName<RealTensorValue>("effective_perm")),
    // Exponential permeability model
    _exponential_permeability_model(getParam<bool>("exponential_permeability_model")),
    _coeff_b(getParam<Real>("coeff_b")),
    // Darcy-Poiseuille permeability model
    _darcy_poiseuille_permeability_model(getParam<bool>("darcy_poiseuille_permeability_model")),
    _wc(getParam<Real>("wc")),
    _perm_exponent(getParam<Real>("perm_exponent"))
{
  _local_elastic_vector.resize(9);
}

void
ADFarmsComputeSmearedCrackingStressGrads::initQpStatefulProperties()
{
  _crack_damage[_qp] = _initial_crack_damage[_qp];
  _eqstrain_local[_qp] = 0.0;
  _kappa[_qp] = 0.0;
  _crack_rotation[_qp] = RankTwoTensor::Identity();
}

void
ADFarmsComputeSmearedCrackingStressGrads::computeQpStress()
{
  // (0) Elastic strain update
  _elastic_strain[_qp] = _elastic_strain_old[_qp] + _strain_increment[_qp];

  // (1) Retrieve material parameters and compute cracking strain ε₀
  const ADReal E = ElasticityTensorTools::getIsotropicYoungsModulus(_elasticity_tensor[_qp]);
  const ADReal tiny = 1e-14;
  const ADReal eps0 = _cracking_stress[_qp] / (E + tiny);

  // (2) Compute Mazars‐type equivalent strain ε̃ and principal directions
  ADRealVectorValue eps_dir;
  computeCrackStrainAndOrientation(eps_dir);
  ADReal p0 = Macaulay(eps_dir(0),false);
  ADReal p1 = Macaulay(eps_dir(1),false);
  ADReal p2 = Macaulay(eps_dir(2),false);
  ADReal eqstrain_local = std::sqrt(p0 * p0 + p1 * p1 + p2 * p2);
  _eqstrain_local[_qp] = eqstrain_local;

  // History with nonlocal regularization (monotonic): kappa = max(kappa_old, nonlocal_eq)
  _kappa[_qp] = std::fmax(_kappa_old[_qp], _eqstrain_nonlocal[_qp]);
  ADReal kappa = _kappa[_qp];

  // Damage law (ensure safe when kappa ~ eps0)
  ADReal omega = _initial_crack_damage[_qp];
  if (kappa > eps0)
  {
    ADReal term = 1.0 - eps0 / (kappa + tiny) * ((1.0 - _paramA) + _paramA * std::exp(_paramB * (eps0 - kappa)));
    omega = std::fmax(term, omega);
  }

  // Irreversibility
  omega = std::fmax(omega, _crack_damage_old[_qp]);
  // Clamp upper bound
  omega = std::fmin(omega, 0.999999);
  _crack_damage[_qp] = omega;

  // update the local elasticity tensor
  updateLocalElasticityTensor();

  // (6) Assign stress - the tangent is automatically computed by AD
  _stress[_qp] = _local_elasticity_tensor * _elastic_strain[_qp];

  // (7) Finite‐strain rotation if needed
  if (_perform_finite_strain_rotations)
  {
    finiteStrainRotation(); // Remove the 'true' argument
    _crack_rotation[_qp] = _rotation_increment[_qp] * _crack_rotation[_qp];
  }

  // Compute effective permeability
  updatePermeabilityForCracking();
}

void
ADFarmsComputeSmearedCrackingStressGrads::computeCrackStrainAndOrientation(
    ADRealVectorValue & strain_in_crack_dir)
{
  // The rotation tensor is ordered such that directions for pre-existing cracks appear first
  // in the list of columns.  For example, if there is one existing crack, its direction is in the
  // first column in the rotation tensor.

  std::vector<ADReal> eigval(3, 0.0);
  ADRankTwoTensor eigvec;

  // Extract regular tensor from AD tensor for eigenvalue calculation
  ADRankTwoTensor elastic_strain_nonad = MetaPhysicL::raw_value(_elastic_strain[_qp]);
  elastic_strain_nonad.symmetricEigenvaluesEigenvectors(eigval, eigvec);

  // If the elastic strain is beyond the cracking strain, save the eigen vectors as
  // the rotation tensor. Reverse their order so that the third principal strain
  // (most tensile) will correspond to the first crack.
  ADRankTwoTensor crack_rotation;
  crack_rotation.fillColumn(0, eigvec.column(2));
  crack_rotation.fillColumn(1, eigvec.column(1));
  crack_rotation.fillColumn(2, eigvec.column(0));
  
  _crack_rotation[_qp] = crack_rotation;

  strain_in_crack_dir(0) = eigval[2];
  strain_in_crack_dir(1) = eigval[1];
  strain_in_crack_dir(2) = eigval[0];
}

void
ADFarmsComputeSmearedCrackingStressGrads::updatePermeabilityForCracking()
{
  // If porous flow coupling is not enabled, return
  if (!_porous_flow_coupling)
    return;

  // Get transformation matrix - extract raw values for rotation
  const RankTwoTensor R = MetaPhysicL::raw_value(_crack_rotation[_qp]);

  // Initialize effective permeability new
  ADRankTwoTensor effective_perm_new;

  // Compute the intrinsic permeability
  ADRankTwoTensor perm_intrinsic = _intrinsic_permeability * ADRankTwoTensor::Identity();

  // Initialize effective permeability new
  // exponential permeability model 
  if (_exponential_permeability_model){
    effective_perm_new = perm_intrinsic * std::exp(_crack_damage[_qp] * _coeff_b);
  }
  // darcy-poiseuille permeability model
  else if (_darcy_poiseuille_permeability_model){
    // Compute crack opening
    // wc is the ultimate crack opening
    ADReal w = _crack_damage[_qp] * _wc; 

    // Compute permeability in the damage zone
    ADRankTwoTensor kf = std::pow(w, 2) / (12.0) * ADRankTwoTensor::Identity();

    // Compute permeability
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

void
ADFarmsComputeSmearedCrackingStressGrads::updateLocalElasticityTensor()
{
  
  // ADRealVectorValue stiffness_ratio_local(1.0, 1.0, 1.0);
  ADRealVectorValue stiffness_ratio_local(1.0 - _crack_damage[_qp], 1.0 - _crack_damage[_qp], 1.0 - _crack_damage[_qp]);
  const ADRankTwoTensor & R = _crack_rotation[_qp];

  const ADReal youngs_modulus =
      ElasticityTensorTools::getIsotropicYoungsModulus(_elasticity_tensor[_qp]);

  const ADReal cracking_stress = _cracking_stress[_qp];

  const ADReal & c0 = stiffness_ratio_local(0);
  const ADReal & c1 = stiffness_ratio_local(1);
  const ADReal & c2 = stiffness_ratio_local(2); 

  const ADReal c01 = c0 * c1;
  const ADReal c02 = c0 * c2;
  const ADReal c12 = c1 * c2;

  const ADReal c01_shear_retention = c01;
  const ADReal c02_shear_retention = c02;
  const ADReal c12_shear_retention = c12;

  _local_elastic_vector[0] = _elasticity_tensor[_qp](0, 0, 0, 0) * c0;
  _local_elastic_vector[1] = _elasticity_tensor[_qp](0, 0, 1, 1) * c01;
  _local_elastic_vector[2] = _elasticity_tensor[_qp](0, 0, 2, 2) * c02;
  _local_elastic_vector[3] = _elasticity_tensor[_qp](1, 1, 1, 1) * c1;
  _local_elastic_vector[4] = _elasticity_tensor[_qp](1, 1, 2, 2) * c12;
  _local_elastic_vector[5] = _elasticity_tensor[_qp](2, 2, 2, 2) * c2;
  _local_elastic_vector[6] = _elasticity_tensor[_qp](1, 2, 1, 2) * c12_shear_retention;
  _local_elastic_vector[7] = _elasticity_tensor[_qp](0, 2, 0, 2) * c02_shear_retention;
  _local_elastic_vector[8] = _elasticity_tensor[_qp](0, 1, 0, 1) * c01_shear_retention;

  // Filling with 9 components is sufficient because these are the only nonzero entries
  // for isotropic or orthotropic materials.
  _local_elasticity_tensor.fillFromInputVector(_local_elastic_vector,
                                              ADRankFourTensor::symmetric9);

  // Rotate the modified elasticity tensor back into global coordinates
  _local_elasticity_tensor.rotate(R);
}

ADReal
ADFarmsComputeSmearedCrackingStressGrads::Macaulay(const ADReal x, const bool deriv)
{
  if (deriv)
    return x > 0 ? 1 : 0;
  return 0.5 * (x + std::abs(x));
}