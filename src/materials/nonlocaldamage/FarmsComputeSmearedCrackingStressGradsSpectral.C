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

#include "FarmsComputeSmearedCrackingStressGradsSpectral.h"
#include "ElasticityTensorTools.h"
#include "StressUpdateBase.h"
#include "Conversion.h"

registerMooseObject("farmsApp", FarmsComputeSmearedCrackingStressGradsSpectral);

InputParameters
FarmsComputeSmearedCrackingStressGradsSpectral::validParams()
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
  // Gradient damage coupling term: 0.5 * h * (e - e_tilde)^2
  params.addParam<Real>("h", 0.0, "Gradient-damage coupling modulus h for 0.5*h*(e - e_tilde)^2");
  params.addParam<Real>("fd_delta", 1e-8, "Finite-difference step for equivalent strain derivatives");
  return params;
}

FarmsComputeSmearedCrackingStressGradsSpectral::FarmsComputeSmearedCrackingStressGradsSpectral(const InputParameters & parameters)
  : ComputeMultipleInelasticStress(parameters),
    // strain energy density outputs (match naming from SmallDeformationIsotropicElasticity)
    _psie(declareProperty<Real>(_base_name + std::string("psie"))),
    _psie_active(declareProperty<Real>(_base_name + std::string("psie_active"))),
    // energies for dissipation tracking
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
    _perm_exponent(getParam<Real>("perm_exponent")),
  _solid_bulk_compliance_damaged(declareProperty<Real>("solid_bulk_compliance_damaged")),
  _h(getParam<Real>("h")),
  _fd_delta(getParam<Real>("fd_delta"))
{
  _local_elastic_vector.resize(9);
}

void
FarmsComputeSmearedCrackingStressGradsSpectral::initQpStatefulProperties()
{
  _crack_damage[_qp] = _initial_crack_damage[_qp];
  _eqstrain_local[_qp] = 0.0;
  _kappa[_qp] = 0.0;
  _crack_rotation[_qp] = RankTwoTensor::Identity();
  // initialize energy accumulators
  _accumulated_elastic_energy[_qp] = 0.0;
  _instant_elastic_energy[_qp] = 0.0;
  _fracture_energy[_qp] = 0.0;
}

void
FarmsComputeSmearedCrackingStressGradsSpectral::computeQpStress()
{
  // (0) Elastic strain update
  _elastic_strain[_qp] = _elastic_strain_old[_qp] + _strain_increment[_qp];

  // Isotropic elastic constants and helpers
  const Real E  = ElasticityTensorTools::getIsotropicYoungsModulus(_elasticity_tensor[_qp]);
  const Real G  = ElasticityTensorTools::getIsotropicShearModulus(_elasticity_tensor[_qp]);
  const Real K  = ElasticityTensorTools::getIsotropicBulkModulus(_elasticity_tensor[_qp]);
  const Real lambda = K - 2.0 * G / 3.0;
  const RankTwoTensor I2(RankTwoTensor::initIdentity);
  const RankFourTensor I4(RankFourTensor::initIdentity);
  const RankFourTensor I4_sym(RankFourTensor::initIdentitySymmetricFour);

  // Cracking strain eps0 and small safeguard
  const Real tiny = 1e-14;
  const Real eps0 = _cracking_stress[_qp] / (E + tiny);

  // (1) Mazars-type equivalent strain (positive principal strains only)
  RealVectorValue prin_eps;
  computeCrackStrainAndOrientation(prin_eps);
  const Real p0 = Macaulay(prin_eps(0), /*deriv=*/false);
  const Real p1 = Macaulay(prin_eps(1), /*deriv=*/false);
  const Real p2 = Macaulay(prin_eps(2), /*deriv=*/false);
  const Real eqstrain_local = std::sqrt(p0 * p0 + p1 * p1 + p2 * p2);
  _eqstrain_local[_qp] = eqstrain_local;

  // (2) Nonlocal history (monotone)
  _kappa[_qp] = std::fmax(_kappa_old[_qp], _eqstrain_nonlocal[_qp]);
  const Real kappa = _kappa[_qp];

  // (3) Damage evolution: omega(kappa) with irreversibility
  Real omega = _initial_crack_damage[_qp];
  if (kappa > eps0)
  {
    const Real term =
        1.0 - eps0 / (kappa + tiny) * ((1.0 - _paramA) + _paramA * std::exp(_paramB * (eps0 - kappa)));
    omega = std::fmax(term, omega);
  }
  omega = std::fmax(omega, _crack_damage_old[_qp]);   // irreversibility
  omega = std::fmin(omega, 0.999999);                 // cap
  _crack_damage[_qp] = omega;

  // Degradation factor
  const Real g = std::max(0.0, 1.0 - omega);

  // (4) Spectral split on strain (stress/strain share eigenvectors under isotropy)
  const Real tr_eps = _elastic_strain[_qp].trace();
  const Real tr_eps_pos = tr_eps > 0.0 ? tr_eps : 0.0;

  // Positive part of strain using spectral decomposition (ε⁺)
  const RankTwoTensor strain_pos = spectralPositivePart(_elastic_strain[_qp]);

  // Intact stress: σ_intact = C : ε
  const RankTwoTensor stress_intact = _elasticity_tensor[_qp] * _elastic_strain[_qp];

  // Positive/negative stress parts:
  //   σ⁺ = λ <tr ε> I + 2G ε⁺,  σ⁻ = σ_intact - σ⁺
  const RankTwoTensor stress_pos = lambda * tr_eps_pos * I2 + 2.0 * G * strain_pos;
  const RankTwoTensor stress_neg = stress_intact - stress_pos;

  // (5) Final stress: degrade tensile part only
  _stress[_qp] = g * stress_pos + stress_neg;

  // (5b) Strain energy density using spectral split, aligned with SmallDeformationIsotropicElasticity
  // Intact energy density: 0.5*lambda*(tr eps)^2 + G * eps:eps
  const Real psie_intact = 0.5 * lambda * tr_eps * tr_eps + G * _elastic_strain[_qp].doubleContraction(_elastic_strain[_qp]);
  // Active part: 0.5*lambda*<tr eps>^2 + G * eps_pos:eps_pos
  const Real psie_active = 0.5 * lambda * tr_eps_pos * tr_eps_pos +
                           G * strain_pos.doubleContraction(strain_pos);
  const Real psie_inactive = psie_intact - psie_active;
  _psie_active[_qp] = psie_active;
  _psie[_qp] = g * psie_active + psie_inactive;

  // (5c) Add gradient-damage energy 0.5 * h * (e - e_tilde)^2
  const Real e_local = _eqstrain_local[_qp];
  const Real e_tilde = _eqstrain_nonlocal[_qp];
  if (_h > 0.0)
    _psie[_qp] += 0.5 * _h * (e_local - e_tilde) * (e_local - e_tilde);

  // (6) Consistent tangent (matches NDSmallDeformationIsotropicElasticity spectral path)
  // C_intact = K I⊗I + 2G (I4_sym − 1/3 I⊗I)
  RankFourTensor C_intact = K * I4 + 2.0 * G * (I4_sym - I4 / 3.0);

  // Positive projector P⁺ from ε (RankTwoTensor API builds it from eigenmodes)
  RankTwoTensor   eigvecs;
  std::vector<Real> eigvals(LIBMESH_DIM);
  RankFourTensor P_pos = _elastic_strain[_qp].positiveProjectionEigenDecomposition(eigvals, eigvecs);

  // Volumetric Heaviside H(tr ε)
  const Real H_tr = (tr_eps > 0.0) ? 1.0 : 0.0;

  // C_pos = λ H(tr ε) I⊗I + 2G P⁺
  const RankFourTensor C_pos = lambda * H_tr * I4 + 2.0 * G * P_pos;

  // Final Jacobian base: C = C_intact + (g - 1) C_pos
  _Jacobian_mult[_qp] = C_intact + (g - 1.0) * C_pos;

  // (6b) Add gradient-damage stress and tangent per Appendix B
  if (_h > 0.0)
  {
    // Derivatives of equivalent strain w.r.t. epsilon (FD-based)
    RankTwoTensor de_dE;          // gradient tensor
    RankFourTensor d2e_dEdE;      // Hessian tensor
    equivalentStrainDerivativesFD(_elastic_strain[_qp], _fd_delta, de_dE, d2e_dEdE);

    // Stress contribution: h (e - e_tilde) de/dE
    const Real diff = e_local - e_tilde;
    _stress[_qp] += _h * diff * de_dE;

    // Tangent contribution: h (de/dE ⊗ de/dE) + h (e - e_tilde) d^2e/dE^2
    _Jacobian_mult[_qp] += _h * de_dE.outerProduct(de_dE) + _h * diff * d2e_dEdE;
  }

  // (7) Finite-strain rotation if requested (keeps your existing update)
  if (_perform_finite_strain_rotations)
  {
    finiteStrainRotation();
    _crack_rotation[_qp] = _rotation_increment[_qp] * _crack_rotation[_qp];
  }

  // (9) Mechanical energy bookkeeping (Chunhui):
  // dEa = 1/2 * (sigma_old : dE + sigma_new : dE)  [trapezoidal work increment]
  // Ei  = psi_e (model free energy density with spectral split)
  // Ea  = Ea_old + dEa
  // Fracture energy Ediss = Ea - Ei = Ea - psi_e
  const Real dEa = 0.5 * (_stress_old[_qp].doubleContraction(_strain_increment[_qp]) +
                          _stress[_qp].doubleContraction(_strain_increment[_qp]));
  const Real Ea = _accumulated_elastic_energy_old[_qp] + dEa;
  _accumulated_elastic_energy[_qp] = Ea;
  // Instantaneous stored elastic energy density consistent with spectral degradation
  _instant_elastic_energy[_qp] = _psie[_qp];
  // Cumulative fracture/dissipated energy density (clamped to avoid tiny negatives)
  _fracture_energy[_qp] = std::max(Ea - _psie[_qp], 0.0);

  // Update solid bulk compliance
  updateSolidBulkCompliance();

  // (8) Permeability update stays the same and uses _crack_damage, _crack_rotation
  updatePermeabilityForCracking();
}

// Equivalent strain from tensor: e = sqrt(sum_i <eps_i>^2)
Real
FarmsComputeSmearedCrackingStressGradsSpectral::equivalentStrainFromTensor(const RankTwoTensor & eps) const
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

// Symmetric basis tensor for (i,j)
RankTwoTensor
FarmsComputeSmearedCrackingStressGradsSpectral::symmetricBasis(unsigned int i, unsigned int j) const
{
  RankTwoTensor B; // initialized to zero
  if (i == j)
    B(i, j) = 1.0;
  else
  {
    B(i, j) = 0.5;
    B(j, i) = 0.5;
  }
  return B;
}

// Finite-difference derivatives of equivalent strain wrt epsilon
void
FarmsComputeSmearedCrackingStressGradsSpectral::equivalentStrainDerivativesFD(const RankTwoTensor & eps,
                                                                              Real delta,
                                                                              RankTwoTensor & grad,
                                                                              RankFourTensor & hess) const
{
  // Initialize outputs
  grad.zero();
  hess.zero();

  // Build list of symmetric index pairs (i<=j)
  std::vector<std::pair<unsigned, unsigned>> idx;
  idx.reserve(6);
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = i; j < 3; ++j)
      idx.emplace_back(i, j);

  // Central differences for gradient
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

  // Second derivatives (mixed) via 4-point stencil
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
      // Accumulate d2 * (Bij ⊗ Bkl)
      hess += d2 * Bij.outerProduct(Bkl);
    }
  }
}

void
FarmsComputeSmearedCrackingStressGradsSpectral::computeCrackStrainAndOrientation(
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
FarmsComputeSmearedCrackingStressGradsSpectral::updateSolidBulkCompliance()
{
  
  // If porous flow coupling is not enabled, return
  if (!_porous_flow_coupling)
    return;

  const Real K  = ElasticityTensorTools::getIsotropicBulkModulus(_elasticity_tensor[_qp]);
  // Update damaged solid bulk compliance C_s(d) = 1 / (g(d) * K)
  // Use a small floor on g to avoid division by zero when damage is nearly complete.
  const Real g_eff = std::max((1-_crack_damage[_qp]), 1e-12);
  // K may be spatially varying; evaluate at current qp
  const Real K_eff = K * g_eff;
  // Declare/update property lazily via reference member
  _solid_bulk_compliance_damaged[_qp] = 1.0 / std::max(K_eff, 1e-24);  

}

void
FarmsComputeSmearedCrackingStressGradsSpectral::updatePermeabilityForCracking()
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

Real
FarmsComputeSmearedCrackingStressGradsSpectral::Macaulay(const Real x, const bool deriv)
{
  if (deriv)
    return x > 0 ? 1 : 0;
  return 0.5 * (x + std::abs(x));
}

// Positive part of a symmetric tensor via spectral decomposition
RankTwoTensor
FarmsComputeSmearedCrackingStressGradsSpectral::spectralPositivePart(const RankTwoTensor & A) const
{
  RankTwoTensor eigvecs;
  std::vector<Real> eigvals(LIBMESH_DIM);
  A.symmetricEigenvaluesEigenvectors(eigvals, eigvecs);

  // Build diagonal of positive eigenvalues
  RankTwoTensor eigvals_pos;
  std::vector<Real> pos(3);
  for (unsigned i = 0; i < 3; ++i)
    pos[i] = eigvals[i] > 0.0 ? eigvals[i] : 0.0;
  eigvals_pos.fillFromInputVector(pos);  // fills diagonal entries

  // Recompose ε⁺ = Q diag(<λ_i>) Q^T
  return eigvecs * eigvals_pos * eigvecs.transpose();
}
