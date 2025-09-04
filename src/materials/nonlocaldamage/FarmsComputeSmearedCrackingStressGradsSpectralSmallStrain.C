//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain.h"
#include "ElasticityTensorTools.h"
#include "SymmetricRankFourTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("farmsApp", FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain);

InputParameters
FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain::validParams()
{
  InputParameters params = ComputeGeneralStressBase::validParams();
  params.addClassDescription("Small-strain spectral smeared-crack model with energy bookkeeping and optional gradient-damage coupling");
  params.addRequiredCoupledVar("cracking_stress",
                               "The stress threshold beyond which cracking occurs. Negative values prevent cracking.");
  params.addRequiredCoupledVar("nonlocal_eqstrain",
                               "The nonlocal equivalent strain used in the damage evolution law");
  params.addRequiredParam<Real>("paramA", "parameter used in the damage evolution law");
  params.addRequiredParam<Real>("paramB", "parameter used in the damage evolution law");
  params.addRequiredCoupledVar("initial_crack_damage",
                               "Initial damage for crack_damage material property");

  params.addParam<bool>("porous_flow_coupling", false, "Enable porous flow coupling");
  params.addParam<Real>("intrinsic_permeability", 5e-19, "Intrinsic permeability in m^2");
  params.addParam<bool>("exponential_permeability_model", false,
                        "Use an exponential function for the effective permeability");
  params.addParam<Real>("coeff_b", -1.0,
                        "Coefficient for the exponential function in the effective permeability");
  params.addParam<bool>("darcy_poiseuille_permeability_model", false,
                        "Use Darcy-Poiseuille model for the effective permeability");
  params.addParam<Real>("wc", -1.0, "ultimate crack width for Darcy-Poiseuille model");
  params.addParam<Real>("perm_exponent", -1.0,
                        "Exponent for the Darcy-Poiseuille model for the effective permeability");
  params.addParam<Real>("h", 0.0, "Gradient-damage coupling modulus h for 0.5*h*(e - e_tilde)^2");
  params.addParam<Real>("fd_delta", 1e-8, "Finite-difference step for equivalent strain derivatives");
  return params;
}

FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain::
    FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain(const InputParameters & parameters)
  : ComputeGeneralStressBase(parameters),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")),
  _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _stress_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "stress")),
    _psie(declareProperty<Real>(_base_name + std::string("psie"))),
    _psie_active(declareProperty<Real>(_base_name + std::string("psie_active"))),
    _accumulated_elastic_energy(declareProperty<Real>(_base_name + std::string("accumulated_elastic_energy"))),
    _accumulated_elastic_energy_old(getMaterialPropertyOld<Real>(_base_name + std::string("accumulated_elastic_energy"))),
    _instant_elastic_energy(declareProperty<Real>(_base_name + std::string("instant_elastic_energy"))),
    _fracture_energy(declareProperty<Real>(_base_name + std::string("fracture_energy"))),
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
    _porous_flow_coupling(getParam<bool>("porous_flow_coupling")),
    _intrinsic_permeability(getParam<Real>("intrinsic_permeability")),
    _effective_perm(declareProperty<RealTensorValue>("effective_perm")),
    _effective_perm_old(getMaterialPropertyOldByName<RealTensorValue>("effective_perm")),
    _exponential_permeability_model(getParam<bool>("exponential_permeability_model")),
    _coeff_b(getParam<Real>("coeff_b")),
    _darcy_poiseuille_permeability_model(getParam<bool>("darcy_poiseuille_permeability_model")),
    _wc(getParam<Real>("wc")),
    _perm_exponent(getParam<Real>("perm_exponent")),
    _solid_bulk_compliance_damaged(declareProperty<Real>("solid_bulk_compliance_damaged")),
    _h(getParam<Real>("h")),
    _fd_delta(getParam<Real>("fd_delta")),
    _strain_increment(declareProperty<RankTwoTensor>("strain_increment"))
{
}

void
FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain::initQpStatefulProperties()
{
  _crack_damage[_qp] = _initial_crack_damage[_qp];
  _eqstrain_local[_qp] = 0.0;
  _kappa[_qp] = 0.0;
  _crack_rotation[_qp] = RankTwoTensor::Identity();
  _accumulated_elastic_energy[_qp] = 0.0;
  _instant_elastic_energy[_qp] = 0.0;
  _fracture_energy[_qp] = 0.0;
}

void
FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain::computeQpStress()
{
  // Small-strain mechanical strain is provided by ComputeSmallStrain as _mechanical_strain
  // For small strain in this model, elastic_strain equals mechanical_strain (no inelastic part)
  _elastic_strain[_qp] = _mechanical_strain[_qp];

  // Isotropic elastic constants and helpers
  const Real E  = ElasticityTensorTools::getIsotropicYoungsModulus(_elasticity_tensor[_qp]);
  const Real G  = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);
  const Real K  = ElasticityTensorTools::getIsotropicBulkModulus(_elasticity_tensor[_qp]);
  const Real lambda = K - 2.0 * G / 3.0;
  const RankTwoTensor I2(RankTwoTensor::initIdentity);
  const RankFourTensor I4(RankFourTensor::initIdentity);
  const RankFourTensor I4_sym(RankFourTensor::initIdentitySymmetricFour);

  // Cracking strain eps0 and safeguard
  const Real tiny = 1e-14;
  const Real eps0 = _cracking_stress[_qp] / (E + tiny);

  // Equivalent strain (Mazars) and orientation
  RealVectorValue prin_eps;
  computeCrackStrainAndOrientation(prin_eps);
  const Real p0 = Macaulay(prin_eps(0), /*deriv=*/false);
  const Real p1 = Macaulay(prin_eps(1), /*deriv=*/false);
  const Real p2 = Macaulay(prin_eps(2), /*deriv=*/false);
  const Real eqstrain_local = std::sqrt(p0 * p0 + p1 * p1 + p2 * p2);
  _eqstrain_local[_qp] = eqstrain_local;

  // Nonlocal history
  _kappa[_qp] = std::fmax(_kappa_old[_qp], _eqstrain_nonlocal[_qp]);
  const Real kappa = _kappa[_qp];

  // Damage evolution with irreversibility
  Real omega = _initial_crack_damage[_qp];
  if (kappa > eps0)
  {
    const Real term = 1.0 - eps0 / (kappa + tiny) * ((1.0 - _paramA) + _paramA * std::exp(_paramB * (eps0 - kappa)));
    omega = std::fmax(term, omega);
  }
  omega = std::fmax(omega, _crack_damage_old[_qp]);
  omega = std::fmin(omega, 0.999999);
  _crack_damage[_qp] = omega;
  const Real g = std::max(0.0, 1.0 - omega);

  // Spectral split on strain
  const Real tr_eps = _elastic_strain[_qp].trace();
  const Real tr_eps_pos = tr_eps > 0.0 ? tr_eps : 0.0;
  const RankTwoTensor strain_pos = spectralPositivePart(_elastic_strain[_qp]);

  // Intact and split stresses
  const RankTwoTensor stress_intact = _elasticity_tensor[_qp] * _elastic_strain[_qp];
  const RankTwoTensor stress_pos = lambda * tr_eps_pos * I2 + 2.0 * G * strain_pos;
  const RankTwoTensor stress_neg = stress_intact - stress_pos;

  // Final stress
  _stress[_qp] = g * stress_pos + stress_neg;

  // Energy density with spectral split
  const Real psie_intact = 0.5 * lambda * tr_eps * tr_eps +
                           G * _elastic_strain[_qp].doubleContraction(_elastic_strain[_qp]);
  const Real psie_active = 0.5 * lambda * tr_eps_pos * tr_eps_pos +
                           G * strain_pos.doubleContraction(strain_pos);
  const Real psie_inactive = psie_intact - psie_active;
  _psie_active[_qp] = psie_active;
  _psie[_qp] = g * psie_active + psie_inactive;

  // Gradient-damage energy
  const Real e_local = _eqstrain_local[_qp];
  const Real e_tilde = _eqstrain_nonlocal[_qp];
  if (_h > 0.0)
    _psie[_qp] += 0.5 * _h * (e_local - e_tilde) * (e_local - e_tilde);

  // Tangent (Jacobian multiplier)
  RankFourTensor C_intact = K * I4 + 2.0 * G * (I4_sym - I4 / 3.0);
  RankTwoTensor   eigvecs;
  std::vector<Real> eigvals(LIBMESH_DIM);
  RankFourTensor P_pos = _elastic_strain[_qp].positiveProjectionEigenDecomposition(eigvals, eigvecs);
  const Real H_tr = (tr_eps > 0.0) ? 1.0 : 0.0;
  const RankFourTensor C_pos = lambda * H_tr * I4 + 2.0 * G * P_pos;
  _Jacobian_mult[_qp] = C_intact + (g - 1.0) * C_pos;

  // Gradient-damage stress and tangent
  if (_h > 0.0)
  {
    RankTwoTensor de_dE;
    RankFourTensor d2e_dEdE;
    equivalentStrainDerivativesFD(_elastic_strain[_qp], _fd_delta, de_dE, d2e_dEdE);
    const Real diff = e_local - e_tilde;
    _stress[_qp] += _h * diff * de_dE;
    _Jacobian_mult[_qp] += _h * de_dE.outerProduct(de_dE) + _h * diff * d2e_dEdE;
  }

  // Energy bookkeeping using mechanical strain increment from framework
  const RankTwoTensor deps = _mechanical_strain[_qp] - _mechanical_strain_old[_qp];
  const Real dEa = 0.5 * (_stress_old[_qp].doubleContraction(deps) +
                          _stress[_qp].doubleContraction(deps));
  const Real Ea = _accumulated_elastic_energy_old[_qp] + dEa;
  _accumulated_elastic_energy[_qp] = Ea;
  _instant_elastic_energy[_qp] = _psie[_qp];
  _fracture_energy[_qp] = std::max(Ea - _psie[_qp], 0.0);

  // Strain increment
  _strain_increment[_qp] = deps;

  // HM updates
  updateSolidBulkCompliance();
  updatePermeabilityForCracking();
}

Real
FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain::Macaulay(const Real x, const bool deriv)
{
  if (deriv)
    return x > 0 ? 1 : 0;
  return 0.5 * (x + std::abs(x));
}

RankTwoTensor
FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain::spectralPositivePart(const RankTwoTensor & A) const
{
  RankTwoTensor eigvecs;
  std::vector<Real> eigvals(LIBMESH_DIM);
  A.symmetricEigenvaluesEigenvectors(eigvals, eigvecs);

  RankTwoTensor eigvals_pos;
  std::vector<Real> pos(3);
  for (unsigned i = 0; i < 3; ++i)
    pos[i] = eigvals[i] > 0.0 ? eigvals[i] : 0.0;
  eigvals_pos.fillFromInputVector(pos);
  return eigvecs * eigvals_pos * eigvecs.transpose();
}

Real
FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain::equivalentStrainFromTensor(const RankTwoTensor & eps) const
{
  RankTwoTensor eigvecs;
  std::vector<Real> eigvals(LIBMESH_DIM);
  eps.symmetricEigenvaluesEigenvectors(eigvals, eigvecs);
  Real acc = 0.0;
  for (unsigned i = 0; i < 3; ++i)
  {
    const Real pi = (eigvals[i] > 0.0) ? eigvals[i] : 0.0;
    acc += pi * pi;
  }
  return std::sqrt(acc);
}

RankTwoTensor
FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain::symmetricBasis(unsigned int i, unsigned int j) const
{
  RankTwoTensor B;
  if (i == j)
    B(i, j) = 1.0;
  else
  {
    B(i, j) = 0.5;
    B(j, i) = 0.5;
  }
  return B;
}

void
FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain::equivalentStrainDerivativesFD(const RankTwoTensor & eps,
                                                                                           Real delta,
                                                                                           RankTwoTensor & grad,
                                                                                           RankFourTensor & hess) const
{
  grad.zero();
  hess.zero();

  std::vector<std::pair<unsigned, unsigned>> idx;
  idx.reserve(6);
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = i; j < 3; ++j)
      idx.emplace_back(i, j);

  for (const auto & ij : idx)
  {
    const unsigned i = ij.first, j = ij.second;
    const RankTwoTensor Bij = symmetricBasis(i, j);
    RankTwoTensor eps_p = eps + delta * Bij;
    RankTwoTensor eps_m = eps - delta * Bij;
    const Real e_p = equivalentStrainFromTensor(eps_p);
    const Real e_m = equivalentStrainFromTensor(eps_m);
    const Real d = 0.5 * (e_p - e_m) / delta;
    grad(i, j) = d;
    grad(j, i) = d;
  }

  for (const auto & ij : idx)
  {
    const unsigned i = ij.first, j = ij.second;
    const RankTwoTensor Bij = symmetricBasis(i, j);
    for (const auto & kl : idx)
    {
      const unsigned k = kl.first, l = kl.second;
      const RankTwoTensor Bkl = symmetricBasis(k, l);

      const RankTwoTensor eps_pp = eps + delta * Bij + delta * Bkl;
      const RankTwoTensor eps_pm = eps + delta * Bij - delta * Bkl;
      const RankTwoTensor eps_mp = eps - delta * Bij + delta * Bkl;
      const RankTwoTensor eps_mm = eps - delta * Bij - delta * Bkl;

      const Real e_pp = equivalentStrainFromTensor(eps_pp);
      const Real e_pm = equivalentStrainFromTensor(eps_pm);
      const Real e_mp = equivalentStrainFromTensor(eps_mp);
      const Real e_mm = equivalentStrainFromTensor(eps_mm);

      const Real d2 = (e_pp - e_pm - e_mp + e_mm) / (4.0 * delta * delta);
      hess += d2 * Bij.outerProduct(Bkl);
    }
  }
}

void
FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain::computeCrackStrainAndOrientation(
    RealVectorValue & strain_in_crack_dir)
{
  std::vector<Real> eigval(3, 0.0);
  RankTwoTensor eigvec;

  _elastic_strain[_qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);

  _crack_rotation[_qp].fillColumn(0, eigvec.column(2));
  _crack_rotation[_qp].fillColumn(1, eigvec.column(1));
  _crack_rotation[_qp].fillColumn(2, eigvec.column(0));

  strain_in_crack_dir(0) = eigval[2];
  strain_in_crack_dir(1) = eigval[1];
  strain_in_crack_dir(2) = eigval[0];
}

void
FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain::updateSolidBulkCompliance()
{
  if (!_porous_flow_coupling)
    return;
  const Real K  = ElasticityTensorTools::getIsotropicBulkModulus(_elasticity_tensor[_qp]);
  const Real g_eff = std::max((1 - _crack_damage[_qp]), 1e-12);
  const Real K_eff = K * g_eff;
  _solid_bulk_compliance_damaged[_qp] = 1.0 / std::max(K_eff, 1e-24);
}

void
FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain::updatePermeabilityForCracking()
{
  if (!_porous_flow_coupling)
    return;

  const RankTwoTensor & R = _crack_rotation[_qp];
  RankTwoTensor effective_perm_new;
  RankTwoTensor perm_intrinsic = _intrinsic_permeability * RankTwoTensor::Identity();

  if (_exponential_permeability_model)
    effective_perm_new = perm_intrinsic * std::exp(_crack_damage[_qp] * _coeff_b);
  else if (_darcy_poiseuille_permeability_model)
  {
    Real w = _crack_damage[_qp] * _wc;
    RankTwoTensor kf = std::pow(w, 2) / 12.0 * RankTwoTensor::Identity();
    effective_perm_new = perm_intrinsic + std::pow(_crack_damage[_qp], _perm_exponent) * (kf - perm_intrinsic);
  }
  else
    mooseError("Unknown permeability model type.");

  effective_perm_new.rotate(R);
  _effective_perm[_qp] = effective_perm_new;
}
