//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "SmallDeformationIsotropicElasticityHM.h"
#include "RaccoonUtils.h"

registerMooseObject("farmsApp", SmallDeformationIsotropicElasticityHM);

InputParameters
SmallDeformationIsotropicElasticityHM::validParams()
{
  InputParameters params = SmallDeformationElasticityModel::validParams();
  params.addClassDescription("Isotropic elasticity under small strain asumptions.");

  params.addRequiredParam<MaterialPropertyName>("bulk_modulus", "The bulk modulus $K$");
  params.addRequiredParam<MaterialPropertyName>("shear_modulus", "The shear modulus $G$");

  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  params.addParam<MaterialPropertyName>(
      "strain_energy_density",
      "psie",
      "Name of the strain energy density computed by this material model");
  params.addParam<MaterialPropertyName>("degradation_function", "g", "The degradation function");
  params.addParam<MooseEnum>(
      "decomposition", MooseEnum("NONE SPECTRAL VOLDEV", "NONE"), "The decomposition method");

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

SmallDeformationIsotropicElasticityHM::SmallDeformationIsotropicElasticityHM(
    const InputParameters & parameters)
  : SmallDeformationElasticityModel(parameters),
    DerivativeMaterialPropertyNameInterface(),
    _K(getADMaterialPropertyByName<Real>(prependBaseName("bulk_modulus", true))),
    _G(getADMaterialPropertyByName<Real>(prependBaseName("shear_modulus", true))),

    _d_name(getVar("phase_field", 0)->name()),

    // The strain energy density and its derivatives
    _psie_name(prependBaseName("strain_energy_density", true)),
    _psie(declareADProperty<Real>(_psie_name)),
    _psie_active(declareADProperty<Real>(_psie_name + "_active")),
    _dpsie_dd(declareADProperty<Real>(derivativePropertyName(_psie_name, {_d_name}))),

    // The degradation function and its derivatives
    _g_name(prependBaseName("degradation_function", true)),
    _g(getADMaterialProperty<Real>(_g_name)),
    _dg_dd(getADMaterialProperty<Real>(derivativePropertyName(_g_name, {_d_name}))),

    _decomposition(getParam<MooseEnum>("decomposition").getEnum<Decomposition>()),

    //porous flow coupling
    _crack_rotation(declareADProperty<RankTwoTensor>("crack_rotation")),
    _crack_rotation_old(getMaterialPropertyOldByName<RankTwoTensor>("crack_rotation")),
    _effective_perm(declareADProperty<RealTensorValue>("effective_perm")),
    _effective_perm_old(getMaterialPropertyOldByName<RealTensorValue>("effective_perm")),
    _porous_flow_coupling(getParam<bool>("porous_flow_coupling")),
    _intrinsic_permeability(getParam<Real>("intrinsic_permeability")),
    
    // Exponential permeability model
    _exponential_permeability_model(getParam<bool>("exponential_permeability_model")),
    _coeff_b(getParam<Real>("coeff_b")),
    // Darcy-Poiseuille permeability model
    _darcy_poiseuille_permeability_model(getParam<bool>("darcy_poiseuille_permeability_model")),
    _wc(getParam<Real>("wc")),
    _perm_exponent(getParam<Real>("perm_exponent")),

    // The phase-field variable
    _d(adCoupledValue("phase_field"))
{
}

ADRankTwoTensor
SmallDeformationIsotropicElasticityHM::computeStress(const ADRankTwoTensor & strain)
{
  ADRankTwoTensor stress;

  if (_decomposition == Decomposition::none)
    stress = computeStressNoDecomposition(strain);
  else if (_decomposition == Decomposition::spectral)
    stress = computeStressSpectralDecomposition(strain);
  else if (_decomposition == Decomposition::voldev)
    stress = computeStressVolDevDecomposition(strain);
  else
    paramError("decomposition", "Unsupported decomposition type.");

  return stress;
}

ADRankTwoTensor
SmallDeformationIsotropicElasticityHM::computeStressNoDecomposition(const ADRankTwoTensor & strain)
{
  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  ADRankTwoTensor stress_intact = _K[_qp] * strain.trace() * I2 + 2 * _G[_qp] * strain.deviatoric();
  ADRankTwoTensor stress = _g[_qp] * stress_intact;

  _psie_active[_qp] = 0.5 * stress_intact.doubleContraction(strain);
  _psie[_qp] = _g[_qp] * _psie_active[_qp];
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  return stress;
}

ADRankTwoTensor
SmallDeformationIsotropicElasticityHM::computeStressSpectralDecomposition(
    const ADRankTwoTensor & strain)
{
  const ADReal lambda = _K[_qp] - 2 * _G[_qp] / LIBMESH_DIM;
  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  ADReal strain_tr = strain.trace();
  ADReal strain_tr_pos = RaccoonUtils::Macaulay(strain_tr);

  // Spectral decomposition
  ADRankTwoTensor strain_pos = RaccoonUtils::spectralDecomposition(strain);

  // Stress
  ADRankTwoTensor stress_intact = _K[_qp] * strain.trace() * I2 + 2 * _G[_qp] * strain.deviatoric();
  ADRankTwoTensor stress_pos = lambda * strain_tr_pos * I2 + 2 * _G[_qp] * strain_pos;
  ADRankTwoTensor stress_neg = stress_intact - stress_pos;
  ADRankTwoTensor stress = _g[_qp] * stress_pos + stress_neg;

  // Strain energy density
  ADReal psie_intact =
      0.5 * lambda * strain_tr * strain_tr + _G[_qp] * strain.doubleContraction(strain);
  _psie_active[_qp] = 0.5 * lambda * strain_tr_pos * strain_tr_pos +
                      _G[_qp] * strain_pos.doubleContraction(strain_pos);
  ADReal psie_inactive = psie_intact - _psie_active[_qp];
  _psie[_qp] = _g[_qp] * _psie_active[_qp] + psie_inactive;
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  //Porous flow coupling
  /* Compute Principal Strains and Rotation Matrix */
  ADRealVectorValue strain_in_crack_dir; //principal strains
  computeCrackStrainAndOrientation(strain_in_crack_dir);

  // Compute effective permeability
  updatePermeabilityForCracking();

  return stress;
}

ADRankTwoTensor
SmallDeformationIsotropicElasticityHM::computeStressVolDevDecomposition(
    const ADRankTwoTensor & strain)
{
  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);

  // Volumetric-deviatoric decomposition
  ADReal strain_tr = strain.trace();
  ADReal strain_tr_pos = RaccoonUtils::Macaulay(strain_tr);
  ADReal strain_tr_neg = strain_tr - strain_tr_pos;
  ADRankTwoTensor strain_dev = strain.deviatoric();

  // Stress
  ADRankTwoTensor stress_intact = _K[_qp] * strain.trace() * I2 + 2 * _G[_qp] * strain.deviatoric();
  ADRankTwoTensor stress_neg = _K[_qp] * strain_tr_neg * I2;
  ADRankTwoTensor stress_pos = stress_intact - stress_neg;
  ADRankTwoTensor stress = _g[_qp] * stress_pos + stress_neg;

  // Strain energy density
  ADReal psie_intact =
      0.5 * _K[_qp] * strain_tr * strain_tr + _G[_qp] * strain_dev.doubleContraction(strain_dev);
  ADReal psie_inactive = 0.5 * _K[_qp] * strain_tr_neg * strain_tr_neg;
  _psie_active[_qp] = psie_intact - psie_inactive;
  _psie[_qp] = _g[_qp] * _psie_active[_qp] + psie_inactive;
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  return stress;
}

void
SmallDeformationIsotropicElasticityHM::computeCrackStrainAndOrientation(
    ADRealVectorValue & strain_in_crack_dir)
{
  // The rotation tensor is ordered such that directions for pre-existing cracks appear first
  // in the list of columns.  For example, if there is one existing crack, its direction is in the
  // first column in the rotation tensor.

  // If porous flow coupling is not enabled, return
  if (!_porous_flow_coupling)
    return;

  std::vector<ADReal> eigval(3, 0.0);
  ADRankTwoTensor eigvec;

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
SmallDeformationIsotropicElasticityHM::updatePermeabilityForCracking()
{

  //Compute the intrinsic permeability
  ADRankTwoTensor perm_intrinsic = _intrinsic_permeability * ADRankTwoTensor::Identity();

  // If porous flow coupling is not enabled, return
  if (!_porous_flow_coupling){
    _effective_perm[_qp] = perm_intrinsic;
    return;
  }

  // Get transformation matrix
  const ADRankTwoTensor & R = _crack_rotation[_qp];

  // Initialize effective permeability new
  ADRankTwoTensor effective_perm_new;

  // Initialize effective permeability new
  // exponential permeability model 
  if (_exponential_permeability_model){
    effective_perm_new = perm_intrinsic * std::exp( _d[_qp] * _coeff_b );
  }
  // darcy-poiseuille permeability model
  else if (_darcy_poiseuille_permeability_model){
    //Compute crack opening
    //wc is the ultimate crack opening
    ADReal w = _d[_qp] * _wc; 

    //Compute permeability in the damage zone
    ADRankTwoTensor kf = std::pow(w, 2) / (12.0) * ADRankTwoTensor::Identity();

    //Compute permeability
    effective_perm_new = perm_intrinsic + std::pow(_d[_qp], _perm_exponent) * (kf - perm_intrinsic);
  }
  else {
    mooseError("Unknown permeability model type.");
  }

  // Rotate back to global frame
  effective_perm_new.rotate(R);

  // Update effective perm
  _effective_perm[_qp] = effective_perm_new;

}
