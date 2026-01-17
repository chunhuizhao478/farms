//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADSmallDeformationIsotropicElasticityPF.h"
#include "RaccoonUtils.h"

registerMooseObject("farmsApp", ADSmallDeformationIsotropicElasticityPF);

InputParameters
ADSmallDeformationIsotropicElasticityPF::validParams()
{
  InputParameters params = SmallDeformationElasticityModel::validParams();
  params.addClassDescription(
      "AD isotropic elasticity with phase-field damage for three-field poroelastodynamics. "
      "Provides stress, strain energy, degradation function, and damage-dependent permeability.");

  params.addRequiredParam<MaterialPropertyName>("bulk_modulus", "The bulk modulus $K$");
  params.addRequiredParam<MaterialPropertyName>("shear_modulus", "The shear modulus $G$");

  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");

  params.addParam<MaterialPropertyName>(
      "strain_energy_density",
      "psie",
      "Name of the strain energy density computed by this material model");
  params.addParam<MaterialPropertyName>("strain_energy_density_active",
                                        "psie_active",
                                        "Name of the active strain energy density");
  params.addParam<MaterialPropertyName>("strain_energy_density_inactive",
                                        "psie_inactive",
                                        "Name of the inactive strain energy density");
  params.addParam<MaterialPropertyName>(
      "strain_energy_density_derivative",
      "dpsie_dd",
      "Name of the strain energy density derivative w/r/t damage");
  params.addParam<MaterialPropertyName>("degradation_function", "g", "The degradation function");
  params.addParam<MaterialPropertyName>("degradation_function_derivative",
                                        "dg_dd",
                                        "Name of the degradation function derivative w/r/t damage");
  params.addParam<MaterialPropertyName>(
      "degradation_function_second_derivative",
      "d2g_dd2",
      "Name of the degradation function second derivative w/r/t damage");
  params.addParam<MooseEnum>(
      "decomposition", MooseEnum("NONE SPECTRAL VOLDEV", "NONE"), "The decomposition method");

  params.addParam<std::string>("model_type", "AT1", "The type of the model: AT1, AT2");
  params.addRequiredParam<Real>("eta", "Parameter in the degradation function (residual stiffness)");

  params.addParam<bool>("porous_flow_coupling", false, "Enable porous flow coupling");
  params.addParam<Real>("intrinsic_permeability", 5e-19, "Intrinsic permeability in m^2");

  // Exponential permeability model
  params.addParam<bool>("exponential_permeability_model",
                        false,
                        "Use an exponential function for the effective permeability");
  params.addParam<Real>("coeff_b",
                        -1.0,
                        "Coefficient for the exponential function in the effective permeability");

  // Darcy-Poiseuille permeability model
  params.addParam<bool>("darcy_poiseuille_permeability_model",
                        false,
                        "Use Darcy-Poiseuille model for the effective permeability");
  params.addParam<Real>("wc", -1.0, "Ultimate crack width for Darcy-Poiseuille model");
  params.addParam<Real>(
      "perm_exponent", -1.0, "Exponent for the Darcy-Poiseuille model for the effective permeability");

  return params;
}

ADSmallDeformationIsotropicElasticityPF::ADSmallDeformationIsotropicElasticityPF(
    const InputParameters & parameters)
  : SmallDeformationElasticityModel(parameters),
    DerivativeMaterialPropertyNameInterface(),
    _K(getADMaterialPropertyByName<Real>(prependBaseName("bulk_modulus", true))),
    _G(getADMaterialPropertyByName<Real>(prependBaseName("shear_modulus", true))),
    _d(adCoupledValue("phase_field")),
    _model_type(getParam<std::string>("model_type")),
    _eta(getParam<Real>("eta")),
    _decomposition(getParam<MooseEnum>("decomposition").getEnum<Decomposition>()),
    _porous_flow_coupling(getParam<bool>("porous_flow_coupling")),
    _intrinsic_permeability(getParam<Real>("intrinsic_permeability")),
    _darcy_poiseuille_permeability_model(getParam<bool>("darcy_poiseuille_permeability_model")),
    _wc(getParam<Real>("wc")),
    _perm_exponent(getParam<Real>("perm_exponent")),
    _exponential_permeability_model(getParam<bool>("exponential_permeability_model")),
    _coeff_b(getParam<Real>("coeff_b")),
    // Strain energy density
    _psie(declareADProperty<Real>(getParam<MaterialPropertyName>("strain_energy_density"))),
    _psie_active(
        declareADProperty<Real>(getParam<MaterialPropertyName>("strain_energy_density_active"))),
    _psie_inactive(
        declareADProperty<Real>(getParam<MaterialPropertyName>("strain_energy_density_inactive"))),
    _dpsie_dd(declareADProperty<Real>(
        getParam<MaterialPropertyName>("strain_energy_density_derivative"))),
    // Degradation function
    _g(declareADProperty<Real>(getParam<MaterialPropertyName>("degradation_function"))),
    _dg_dd(declareADProperty<Real>(
        getParam<MaterialPropertyName>("degradation_function_derivative"))),
    _d2g_dd2(declareADProperty<Real>(
        getParam<MaterialPropertyName>("degradation_function_second_derivative"))),
    // Porous flow coupling
    _crack_rotation(declareADProperty<RankTwoTensor>("crack_rotation")),
    _crack_rotation_old(getMaterialPropertyOldByName<RankTwoTensor>("crack_rotation")),
    _effective_perm(declareADProperty<RankTwoTensor>("effective_perm")),
    _effective_perm_old(getMaterialPropertyOldByName<RankTwoTensor>("effective_perm"))
{
  // Parameter validation
  if (_porous_flow_coupling && !_exponential_permeability_model && !_darcy_poiseuille_permeability_model)
    paramError("porous_flow_coupling",
               "Porous flow coupling is enabled, but no permeability model is selected. "
               "Please enable either exponential or Darcy-Poiseuille permeability model.");
  if (_darcy_poiseuille_permeability_model && (_wc <= 0.0 || _perm_exponent <= 0.0))
    paramError("darcy_poiseuille_permeability_model",
               "Darcy-Poiseuille permeability model is enabled, but wc and perm_exponent "
               "must be positive values.");
  if (_exponential_permeability_model && _coeff_b <= 0.0)
    paramError("exponential_permeability_model",
               "Exponential permeability model is enabled, but coeff_b must be a positive value.");
}

void
ADSmallDeformationIsotropicElasticityPF::initQpStatefulProperties()
{
  SmallDeformationElasticityModel::initQpStatefulProperties();

  _psie[_qp] = 0.0;
  _psie_active[_qp] = 0.0;
  _psie_inactive[_qp] = 0.0;
  _dpsie_dd[_qp] = 0.0;
  _g[_qp] = 1.0;
  _dg_dd[_qp] = 0.0;
  _d2g_dd2[_qp] = 0.0;

  // Initialize permeability
  _crack_rotation[_qp] = ADRankTwoTensor::Identity();
  _effective_perm[_qp] = _intrinsic_permeability * ADRankTwoTensor::Identity();
}

ADRankTwoTensor
ADSmallDeformationIsotropicElasticityPF::computeStress(const ADRankTwoTensor & strain)
{
  ADRankTwoTensor stress;

  // Compute degradation function and derivatives
  computeGDerivatives();

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
ADSmallDeformationIsotropicElasticityPF::computeStressNoDecomposition(const ADRankTwoTensor & strain)
{
  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  ADRankTwoTensor stress_intact = _K[_qp] * strain.trace() * I2 + 2 * _G[_qp] * strain.deviatoric();
  ADRankTwoTensor stress = _g[_qp] * stress_intact;

  _psie_active[_qp] = 0.5 * stress_intact.doubleContraction(strain);
  _psie_inactive[_qp] = 0.0;
  _psie[_qp] = _g[_qp] * _psie_active[_qp];
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  return stress;
}

ADRankTwoTensor
ADSmallDeformationIsotropicElasticityPF::computeStressSpectralDecomposition(
    const ADRankTwoTensor & strain)
{
  const ADReal lambda = _K[_qp] - 2 * _G[_qp] / LIBMESH_DIM;
  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  ADReal strain_tr = strain.trace();
  ADReal strain_tr_pos = Macaulay(strain_tr);

  // Spectral decomposition
  ADRankTwoTensor strain_pos = spectralDecomposition(strain);

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
  _psie_inactive[_qp] = psie_inactive;
  _psie[_qp] = _g[_qp] * _psie_active[_qp] + psie_inactive;
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  // Porous flow coupling - compute permeability
  if (_porous_flow_coupling)
  {
    RealVectorValue strain_in_crack_dir;
    computeCrackStrainAndOrientation(strain_in_crack_dir);
    updatePermeabilityForCracking();
  }

  return stress;
}

ADRankTwoTensor
ADSmallDeformationIsotropicElasticityPF::computeStressVolDevDecomposition(
    const ADRankTwoTensor & strain)
{
  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);

  // Volumetric-deviatoric decomposition
  ADReal strain_tr = strain.trace();
  ADReal strain_tr_pos = Macaulay(strain_tr);
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
  _psie_inactive[_qp] = psie_inactive;
  _psie[_qp] = _g[_qp] * _psie_active[_qp] + psie_inactive;
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  return stress;
}

void
ADSmallDeformationIsotropicElasticityPF::computeGDerivatives()
{
  const ADReal d = _d[_qp];

  if (_model_type == "AT2")
  {
    _g[_qp] = std::pow((1 - d), 2) * (1 - _eta) + _eta;
    _dg_dd[_qp] = -2 * (1 - _eta) * (1 - d);
    _d2g_dd2[_qp] = 2 * (1 - _eta);
  }
  else if (_model_type == "AT1")
  {
    _g[_qp] = std::pow((1 - d), 2) * (1 - _eta) + _eta;
    _dg_dd[_qp] = -2 * (1 - _eta) * (1 - d);
    _d2g_dd2[_qp] = 2 * (1 - _eta);
  }
  else
    mooseError("Unknown model type: " + _model_type);
}

void
ADSmallDeformationIsotropicElasticityPF::computeCrackStrainAndOrientation(
    RealVectorValue & strain_in_crack_dir)
{
  if (!_porous_flow_coupling)
    return;

  std::vector<ADReal> eigval(3, 0.0);
  ADRankTwoTensor eigvec;

  _elastic_strain[_qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);

  // Reverse order so that the third principal strain (most tensile) corresponds to the first crack
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
    {
      _crack_rotation[_qp](i, 0) = eigvec(i, 2);
      _crack_rotation[_qp](i, 1) = eigvec(i, 1);
      _crack_rotation[_qp](i, 2) = eigvec(i, 0);
    }

  strain_in_crack_dir(0) = MetaPhysicL::raw_value(eigval[2]);
  strain_in_crack_dir(1) = MetaPhysicL::raw_value(eigval[1]);
  strain_in_crack_dir(2) = MetaPhysicL::raw_value(eigval[0]);
}

void
ADSmallDeformationIsotropicElasticityPF::updatePermeabilityForCracking()
{
  if (!_porous_flow_coupling)
    return;

  // Get transformation matrix
  const ADRankTwoTensor & R = _crack_rotation[_qp];

  // Initialize effective permeability new
  ADRankTwoTensor effective_perm_new;

  // Compute the intrinsic permeability
  ADRankTwoTensor perm_intrinsic = _intrinsic_permeability * ADRankTwoTensor::Identity();

  // Get damage value
  const ADReal d = _d[_qp];

  if (_exponential_permeability_model)
  {
    effective_perm_new = perm_intrinsic * std::exp(d * _coeff_b);
  }
  else if (_darcy_poiseuille_permeability_model)
  {
    // Compute crack opening: wc is the ultimate crack opening
    ADReal w = d * _wc;

    // Compute permeability in the damage zone
    ADRankTwoTensor kf = std::pow(w, 2) / 12.0 * ADRankTwoTensor::Identity();

    // Compute permeability
    effective_perm_new = perm_intrinsic + std::pow(d, _perm_exponent) * (kf - perm_intrinsic);
  }
  else
  {
    mooseError("Unknown permeability model type.");
  }

  // Rotate back to global frame
  effective_perm_new.rotate(R);

  // Update effective perm
  _effective_perm[_qp] = effective_perm_new;
}

ADReal
ADSmallDeformationIsotropicElasticityPF::Macaulay(const ADReal & x)
{
  return 0.5 * (x + std::abs(x));
}

std::vector<ADReal>
ADSmallDeformationIsotropicElasticityPF::Macaulay(const std::vector<ADReal> & v)
{
  std::vector<ADReal> m = v;
  for (auto & x : m)
    x = Macaulay(x);
  return m;
}

ADRankTwoTensor
ADSmallDeformationIsotropicElasticityPF::spectralDecomposition(const ADRankTwoTensor & r2t)
{
  ADRankTwoTensor eigvecs;
  std::vector<ADReal> eigvals(LIBMESH_DIM);
  r2t.symmetricEigenvaluesEigenvectors(eigvals, eigvecs);

  ADRankTwoTensor eigvals_pos;
  eigvals_pos.fillFromInputVector(Macaulay(eigvals));
  return eigvecs * eigvals_pos * eigvecs.transpose();
}
