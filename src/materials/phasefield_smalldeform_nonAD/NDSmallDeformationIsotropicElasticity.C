//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "NDSmallDeformationIsotropicElasticity.h"
#include "RaccoonUtils.h"
#include "SymmetricRankFourTensor.h"

registerMooseObject("farmsApp", NDSmallDeformationIsotropicElasticity);

InputParameters
NDSmallDeformationIsotropicElasticity::validParams()
{
  InputParameters params = NDSmallDeformationElasticityModel::validParams();
  params.addClassDescription("Isotropic elasticity under small strain asumptions.");

  params.addRequiredParam<MaterialPropertyName>("bulk_modulus", "The bulk modulus $K$");
  params.addRequiredParam<MaterialPropertyName>("shear_modulus", "The shear modulus $G$");

  //material property names
  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  params.addParam<MaterialPropertyName>(
      "strain_energy_density",
      "psie",
      "Name of the strain energy density computed by this material model");
  params.addParam<MaterialPropertyName>("strain_energy_density_active",
                                        "psie_active",
                                        "Name of the active strain energy density");
  params.addParam<MaterialPropertyName>(
      "strain_energy_density_derivative",
      "dpsie_dd",
      "Name of the strain energy density derivative w/r/t damage");
  params.addParam<MaterialPropertyName>("degradation_function", "g", "The degradation function");
  params.addParam<MaterialPropertyName>(
      "degradation_function_derivative",
      "dg_dd",
      "Name of the degradation function derivative w/r/t damage");
  params.addParam<MaterialPropertyName>(
      "degradation_function_second_derivative",
      "d2g_dd2",
      "Name of the degradation function second derivative w/r/t damage");
  params.addParam<MooseEnum>(
      "decomposition", MooseEnum("NONE SPECTRAL VOLDEV", "NONE"), "The decomposition method");

  params.addParam<std::string>(
      "model_type",
      "AT1",
      "The type of the model: AT1, AT2, PF_CZM");
  
  params.addRequiredParam<Real>("eta", "Parameter in the degradation function");

  return params;
}

NDSmallDeformationIsotropicElasticity::NDSmallDeformationIsotropicElasticity(
    const InputParameters & parameters)
  : NDSmallDeformationElasticityModel(parameters),
    DerivativeMaterialPropertyNameInterface(),
    _K(getMaterialPropertyByName<Real>(prependBaseName("bulk_modulus", true))),
    _G(getMaterialPropertyByName<Real>(prependBaseName("shear_modulus", true))),

    _d(coupledValue("phase_field")),

    // The strain energy density and its derivatives
    _psie(declareProperty<Real>(getParam<MaterialPropertyName>("strain_energy_density"))),
    _psie_active(declareProperty<Real>(getParam<MaterialPropertyName>(
        "strain_energy_density_active"))),
    _dpsie_dd(declareProperty<Real>(getParam<MaterialPropertyName>(
        "strain_energy_density_derivative"))),

    // The degradation function and its derivatives
    _g(declareProperty<Real>(getParam<MaterialPropertyName>("degradation_function"))),
    _dg_dd(declareProperty<Real>(getParam<MaterialPropertyName>(
        "degradation_function_derivative"))),
    _d2g_dd2(declareProperty<Real>(getParam<MaterialPropertyName>(
        "degradation_function_second_derivative"))),

    // model type
    _model_type(getParam<std::string>("model_type")),

    // Constants
    _eta(getParam<Real>("eta")),

    _decomposition(getParam<MooseEnum>("decomposition").getEnum<Decomposition>())
{
}

RankTwoTensor
NDSmallDeformationIsotropicElasticity::computeStress(const RankTwoTensor & strain)
{
  RankTwoTensor stress;

  // Evaluate g and derivatives
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

RankFourTensor
NDSmallDeformationIsotropicElasticity::computeJacobian(const RankTwoTensor & strain)
{
  RankFourTensor Jacobian;

  // Evaluate g and derivatives
  computeGDerivatives();

  if (_decomposition == Decomposition::none)
    Jacobian = computeJacobianNoDecomposition(strain);
  else if (_decomposition == Decomposition::spectral)
    Jacobian = computeJacobianSpectralDecomposition(strain);
  else if (_decomposition == Decomposition::voldev)
    Jacobian = computeJacobianVolDevDecomposition(strain);
  else
    paramError("decomposition", "Unsupported decomposition type.");

  return Jacobian;
}

RankTwoTensor
NDSmallDeformationIsotropicElasticity::computeStressNoDecomposition(const RankTwoTensor & strain)
{
  const RankTwoTensor I2(RankTwoTensor::initIdentity);
  RankTwoTensor stress_intact = _K[_qp] * strain.trace() * I2 + 2 * _G[_qp] * strain.deviatoric();
  RankTwoTensor stress = _g[_qp] * stress_intact;

  _psie_active[_qp] = 0.5 * stress_intact.doubleContraction(strain);
  _psie[_qp] = _g[_qp] * _psie_active[_qp];
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  return stress;
}

RankTwoTensor
NDSmallDeformationIsotropicElasticity::computeStressSpectralDecomposition(
    const RankTwoTensor & strain)
{
  const Real lambda = _K[_qp] - 2 * _G[_qp] / LIBMESH_DIM;
  const RankTwoTensor I2(RankTwoTensor::initIdentity);
  Real strain_tr = strain.trace();
  Real strain_tr_pos = NDSmallDeformationIsotropicElasticity::Macaulay(strain_tr);

  // Spectral decomposition
  RankTwoTensor strain_pos = NDSmallDeformationIsotropicElasticity::spectralDecomposition(strain);

  // Stress
  RankTwoTensor stress_intact = _K[_qp] * strain.trace() * I2 + 2 * _G[_qp] * strain.deviatoric();
  RankTwoTensor stress_pos = lambda * strain_tr_pos * I2 + 2 * _G[_qp] * strain_pos;
  RankTwoTensor stress_neg = stress_intact - stress_pos;
  RankTwoTensor stress = _g[_qp] * stress_pos + stress_neg;

  // Strain energy density
  Real psie_intact =
      0.5 * lambda * strain_tr * strain_tr + _G[_qp] * strain.doubleContraction(strain);
  _psie_active[_qp] = 0.5 * lambda * strain_tr_pos * strain_tr_pos +
                      _G[_qp] * strain_pos.doubleContraction(strain_pos);
  Real psie_inactive = psie_intact - _psie_active[_qp];
  _psie[_qp] = _g[_qp] * _psie_active[_qp] + psie_inactive;
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  return stress;
}

RankTwoTensor
NDSmallDeformationIsotropicElasticity::computeStressVolDevDecomposition(
    const RankTwoTensor & strain)
{
  const RankTwoTensor I2(RankTwoTensor::initIdentity);

  // Volumetric-deviatoric decomposition
  Real strain_tr = strain.trace();
  Real strain_tr_pos = NDSmallDeformationIsotropicElasticity::Macaulay(strain_tr);
  Real strain_tr_neg = strain_tr - strain_tr_pos;
  RankTwoTensor strain_dev = strain.deviatoric();

  // Stress
  RankTwoTensor stress_intact = _K[_qp] * strain.trace() * I2 + 2 * _G[_qp] * strain.deviatoric();
  RankTwoTensor stress_neg = _K[_qp] * strain_tr_neg * I2;
  RankTwoTensor stress_pos = stress_intact - stress_neg;
  RankTwoTensor stress = _g[_qp] * stress_pos + stress_neg;

  // Strain energy density
  Real psie_intact =
      0.5 * _K[_qp] * strain_tr * strain_tr + _G[_qp] * strain_dev.doubleContraction(strain_dev);
  Real psie_inactive = 0.5 * _K[_qp] * strain_tr_neg * strain_tr_neg;
  _psie_active[_qp] = psie_intact - psie_inactive;
  _psie[_qp] = _g[_qp] * _psie_active[_qp] + psie_inactive;
  _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];

  return stress;
}

//Jacobian of the stress w/r/t strain with decomposition methods
// no decomposition: σ = g·(K tr ε I + 2G ε_dev)
//
// ⇒ J_ijkl = g · [ K δ_ij δ_kl + 2G ( ½(δ_ik δ_jl + δ_il δ_jk) − 1/3 δ_ij δ_kl ) ]
RankFourTensor
NDSmallDeformationIsotropicElasticity::computeJacobianNoDecomposition(
    const RankTwoTensor & strain)
{
  const RankTwoTensor I2(RankTwoTensor::initIdentity);
  RankFourTensor I4 = RankFourTensor(RankFourTensor::initIdentity);
  RankFourTensor I4_sym = RankFourTensor(RankFourTensor::initIdentitySymmetricFour);

  RankFourTensor Jacobian_intact = _K[_qp] * I4 + 2 * _G[_qp] * (I4_sym - I4 / 3.0);
  RankFourTensor Jacobian = _g[_qp] * Jacobian_intact;

  return Jacobian;
}

RankFourTensor
NDSmallDeformationIsotropicElasticity::computeJacobianSpectralDecomposition(
    const RankTwoTensor & strain)
{
  //--------------------------------------------------------------------
  // 1.  Some handy constants and fourth–order identity tensors
  //--------------------------------------------------------------------
  const Real  lambda = _K[_qp] - 2.0 * _G[_qp] / LIBMESH_DIM;

  const RankFourTensor I4 = RankFourTensor(RankFourTensor::initIdentity);
  const RankFourTensor I4_sym = RankFourTensor(RankFourTensor::initIdentitySymmetricFour);

  //--------------------------------------------------------------------
  // 2.  Intact (undegraded) isotropic‑elastic tangent
  //--------------------------------------------------------------------
  RankFourTensor C_intact =
      _K[_qp] * I4 + 2.0 * _G[_qp] * (I4_sym - I4 / 3.0);

  //--------------------------------------------------------------------
  // 3.  Positive‑part projector  P⁺  and volumetric Heaviside term
  //--------------------------------------------------------------------
  RankTwoTensor   eigvecs;
  std::vector<Real> eigvals(LIBMESH_DIM);
  RankFourTensor P_pos =
      strain.positiveProjectionEigenDecomposition(eigvals, eigvecs); // P⁺₍ᵢⱼₖₗ₎

  const Real H_tr = (strain.trace() > 0.0) ? 1.0 : 0.0;              // H(tr ε)

  RankFourTensor C_pos =
      lambda * H_tr * I4           // volumetric part
    + 2.0 * _G[_qp] * P_pos;       // deviatoric spectral part

  //--------------------------------------------------------------------
  // 4.  Final tangent  C = C_intact + (g-1) C_pos
  //--------------------------------------------------------------------
  RankFourTensor Jacobian =
      C_intact + (_g[_qp] - 1.0) * C_pos;

  return Jacobian;
}

RankFourTensor
NDSmallDeformationIsotropicElasticity::computeJacobianVolDevDecomposition(
    const RankTwoTensor & strain)
{
  RankFourTensor Jacobian;
  return Jacobian;
}

//helper function, grab from RaccoonUtils, make it non-AD
Real
NDSmallDeformationIsotropicElasticity::Macaulay(const Real x, const bool deriv)
{
  if (deriv)
    return x > 0 ? 1 : 0;
  return 0.5 * (x + std::abs(x));
}

std::vector<Real>
NDSmallDeformationIsotropicElasticity::Macaulay(const std::vector<Real> & v, const bool deriv)
{
  std::vector<Real> m = v;
  for (auto & x : m)
    x = Macaulay(x, deriv);
  return m;
}

RankTwoTensor
NDSmallDeformationIsotropicElasticity::spectralDecomposition(const RankTwoTensor & r2t)
{
  RankTwoTensor eigvecs;
  std::vector<Real> eigvals(LIBMESH_DIM);
  r2t.symmetricEigenvaluesEigenvectors(eigvals, eigvecs);

  RankTwoTensor eigvals_pos;
  eigvals_pos.fillFromInputVector(Macaulay(eigvals));
  return eigvecs * eigvals_pos * eigvecs.transpose();
}

void
NDSmallDeformationIsotropicElasticity::computeGDerivatives()
{
  if (_model_type == "AT2"){
    _g[_qp] = std::pow( (1 - _d[_qp]) , 2 ) * (1 - _eta) + _eta;
    _dg_dd[_qp] = -2 * (1 - _eta) * (1 - _d[_qp]);
    _d2g_dd2[_qp] = 2 * (1 - _eta);
  }
  else if (_model_type == "AT1"){
    _g[_qp] = std::pow( (1 - _d[_qp]) , 2 ) * (1 - _eta) + _eta;
    _dg_dd[_qp] = -2 * (1 - _eta) * (1 - _d[_qp]);
    _d2g_dd2[_qp] = 2 * (1 - _eta);
  }
  else if (_model_type == "PF_CZM"){
    mooseError("PF_CZM model is not implemented yet.");
  }
  else
    mooseError("Unknown model type: " + _model_type);
}
